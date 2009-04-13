/** @file mul.cpp
 *
 *  Implementation of GiNaC's products of expressions. */

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

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <sstream>

#include "mul.h"
#include "add.h"
#include "operators.h"
#include "matrix.h"
#include "indexed.h"
#include "lst.h"
#include "archive.h"
#include "utils.h"
#include "symbol.h"
#include "compiler.h"
#include "constant.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(mul, expairseq,
  print_func<print_context>(&mul::do_print).
  print_func<print_latex>(&mul::do_print_latex).
  print_func<print_csrc>(&mul::do_print_csrc).
  print_func<print_tree>(&mul::do_print_tree).
  print_func<print_python_repr>(&mul::do_print_python_repr))


//////////
// default constructor
//////////

mul::mul()
{
	tinfo_key = &mul::tinfo_static;
}

//////////
// other constructors
//////////

// public

mul::mul(const ex & lh, const ex & rh)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = _ex1;
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

mul::mul(const exvector & v)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = _ex1;
	construct_from_exvector(v);
	GINAC_ASSERT(is_canonical());
}

mul::mul(const epvector & v)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = _ex1;
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

mul::mul(const epvector & v, const ex & oc, bool do_index_renaming)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = oc;
	construct_from_epvector(v, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

mul::mul(std::auto_ptr<epvector> vp, const ex & oc, bool do_index_renaming)
{
	tinfo_key = &mul::tinfo_static;
	GINAC_ASSERT(vp.get()!=0);
	overall_coeff = oc;
	construct_from_epvector(*vp, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

mul::mul(const ex & lh, const ex & mh, const ex & rh)
{
	tinfo_key = &mul::tinfo_static;
	exvector factors;
	factors.reserve(3);
	factors.push_back(lh);
	factors.push_back(mh);
	factors.push_back(rh);
	overall_coeff = _ex1;
	construct_from_exvector(factors);
	GINAC_ASSERT(is_canonical());
}

//////////
// archiving
//////////

DEFAULT_ARCHIVING(mul)

//////////
// functions overriding virtual functions from base classes
//////////

void mul::print_overall_coeff(const print_context & c,
		const char *mul_sym, bool latex) const
{
	const numeric &coeff = ex_to<numeric>(overall_coeff);
	std::stringstream tstream;
	print_context* tcontext_p;
	if (latex)
		tcontext_p = new print_latex(tstream, c.options);
	else
		tcontext_p = new print_dflt(tstream, c.options);
	//print_context tcontext(tstream, c.options);
	coeff.print(*tcontext_p, 0);
	std::string coeffstr = tstream.str();
	delete tcontext_p;
	if (coeffstr[0] == '-') {
		c.s << '-';
		coeffstr = coeffstr.erase(0, 1);
	}

	bool paranthesis = ((coeffstr.find(' ') != std::string::npos)||
		(coeffstr.find('+') != std::string::npos) ||
		(coeffstr.find('-') != std::string::npos));// ||
		//(coeffstr.find('/') != std::string::npos) ||
		//(coeffstr.find('*') != std::string::npos) ||
		//(coeffstr.find('^') != std::string::npos));
	if (!coeff.is_equal(*_num1_p) &&
		(!coeff.is_equal(*_num_1_p) || coeff.is_parent_pos_char())) {
		if (paranthesis) {
			if (latex)
				c.s << "\\left(";
			else
				c.s << '(';
		}
		c.s<<coeffstr;
		if (paranthesis) {
			if (latex)
				c.s << "\\right)";
			else
				c.s << ')';
		}
		if (latex) {
			c.s << " \\, ";
		} else {
			c.s << mul_sym;
		}
	}
}

/*
void mul::do_print(const print_context & c, unsigned level) const
{
	if (precedence() <= level)
		c.s << '(';

	print_overall_coeff(c, "*");

	epvector::const_iterator it = seq.begin(), itend = seq.end();
	bool first = true;
	while (it != itend) {
		if (!first)
			c.s << '*';
		else
			first = false;
		recombine_pair_to_ex(*it).print(c, precedence());
		++it;
	}

	if (precedence() <= level)
		c.s << ')';
}
*/

void mul::do_print(const print_context & c, unsigned level) const
{
	do_print_rat_func(c, level, false);
}

void mul::do_print_latex(const print_latex & c, unsigned level) const
{
	do_print_rat_func(c, level, true);
}

void mul::print_exvector(const exvector & v, const print_context & c,
		char* sep) const
{
	bool first = true;
	exvector::const_iterator vit = v.begin(), vitend = v.end();
	while (vit != vitend) {
		if (!first)
			c.s << sep;
		else
			first = false;
		vit->print(c, precedence());
		++vit;
	}
}

void mul::do_print_rat_func(const print_context & c, unsigned level,
		bool latex_tags) const
{
	char *sep;
	if (latex_tags) {
		sep = " ";
	} else {
		sep = "*";
	}

	// Separate factors into those with negative numeric exponent
	// and all others
	epvector::const_iterator it = seq.begin(), itend = seq.end();
	exvector neg_powers, others;
	while (it != itend) {
		GINAC_ASSERT(is_exactly_a<numeric>(it->coeff));
		if (ex_to<numeric>(it->coeff).is_real() &&
				ex_to<numeric>(it->coeff).is_negative())
			neg_powers.push_back(recombine_pair_to_ex(expair(it->rest, -(it->coeff))));
		else
			others.push_back(recombine_pair_to_ex(*it));
		++it;
	}

	if (!neg_powers.empty()) {
		// Factors with negative exponent are printed as a fraction
		print_overall_coeff(c, others.size() == 0 ? "" : sep,
				latex_tags);
		if (latex_tags) {
			c.s << "\\frac{";
			if (others.empty())
				c.s<<"1";
			else
				print_exvector(others, c, sep);
			c.s << "}{";
			print_exvector(neg_powers, c, sep);
			c.s << "}";
		} else {
			if (others.empty() && ( overall_coeff.is_equal(_ex1) ||
						overall_coeff.is_equal(_ex_1)))
				c.s<<"1";
			else
				print_exvector(others, c, sep);
			c.s << "/";

			if (neg_powers.size() > 1) {
				c.s<<"(";
			}
			print_exvector(neg_powers, c, sep);
			if (neg_powers.size() > 1) {
				c.s << ")";
			}
		}

	} else {

		print_overall_coeff(c, sep, latex_tags);

		print_exvector(others, c, sep);
	}

}

void mul::do_print_csrc(const print_csrc & c, unsigned level) const
{
	if (precedence() <= level)
		c.s << "(";

	if (!overall_coeff.is_equal(_ex1)) {
		if (overall_coeff.is_equal(_ex_1))
			c.s << "-";
		else {
			overall_coeff.print(c, precedence());
			c.s << "*";
		}
	}

	// Print arguments, separated by "*" or "/"
	epvector::const_iterator it = seq.begin(), itend = seq.end();
	while (it != itend) {

		// If the first argument is a negative integer power, it gets printed as "1.0/<expr>"
		bool needclosingparenthesis = false;
		if (it == seq.begin() && it->coeff.info(info_flags::negint)) {
			if (is_a<print_csrc_cl_N>(c)) {
				c.s << "recip(";
				needclosingparenthesis = true;
			} else
				c.s << "1.0/";
		}

		// If the exponent is 1 or -1, it is left out
		if (it->coeff.is_equal(_ex1) || it->coeff.is_equal(_ex_1))
			it->rest.print(c, precedence());
		else if (it->coeff.info(info_flags::negint))
			// Outer parens around ex needed for broken GCC parser:
			(ex(power(it->rest, -ex_to<numeric>(it->coeff)))).print(c, level);
		else
			// Outer parens around ex needed for broken GCC parser:
			(ex(power(it->rest, ex_to<numeric>(it->coeff)))).print(c, level);

		if (needclosingparenthesis)
			c.s << ")";

		// Separator is "/" for negative integer powers, "*" otherwise
		++it;
		if (it != itend) {
			if (it->coeff.info(info_flags::negint))
				c.s << "/";
			else
				c.s << "*";
		}
	}

	if (precedence() <= level)
		c.s << ")";
}

void mul::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(';
	op(0).print(c);
	for (size_t i=1; i<nops(); ++i) {
		c.s << ',';
		op(i).print(c);
	}
	c.s << ')';
}

bool mul::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::polynomial:
		case info_flags::integer_polynomial:
		case info_flags::cinteger_polynomial:
		case info_flags::rational_polynomial:
		case info_flags::crational_polynomial:
		case info_flags::rational_function: {
			epvector::const_iterator i = seq.begin(), end = seq.end();
			while (i != end) {
				if (!(recombine_pair_to_ex(*i).info(inf)))
					return false;
				++i;
			}
			return overall_coeff.info(inf);
		}
		case info_flags::algebraic: {
			epvector::const_iterator i = seq.begin(), end = seq.end();
			while (i != end) {
				if ((recombine_pair_to_ex(*i).info(inf)))
					return true;
				++i;
			}
			return false;
		}
	}
	return inherited::info(inf);
}

int mul::degree(const ex & s) const
{
	// Sum up degrees of factors
	int deg_sum = 0;
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		if (ex_to<numeric>(i->coeff).is_integer())
			deg_sum += recombine_pair_to_ex(*i).degree(s);
		else {
			if (i->rest.has(s))
				throw std::runtime_error("mul::degree() undefined degree because of non-integer exponent");
		}
		++i;
	}
	return deg_sum;
}

int mul::ldegree(const ex & s) const
{
	// Sum up degrees of factors
	int deg_sum = 0;
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		if (ex_to<numeric>(i->coeff).is_integer())
			deg_sum += recombine_pair_to_ex(*i).ldegree(s);
		else {
			if (i->rest.has(s))
				throw std::runtime_error("mul::ldegree() undefined degree because of non-integer exponent");
		}
		++i;
	}
	return deg_sum;
}

ex mul::coeff(const ex & s, int n) const
{
	exvector coeffseq;
	coeffseq.reserve(seq.size()+1);
	
	if (n==0) {
		// product of individual coeffs
		// if a non-zero power of s is found, the resulting product will be 0
		epvector::const_iterator i = seq.begin(), end = seq.end();
		while (i != end) {
			coeffseq.push_back(recombine_pair_to_ex(*i).coeff(s,n));
			++i;
		}
		coeffseq.push_back(overall_coeff);
		return (new mul(coeffseq))->setflag(status_flags::dynallocated);
	}
	
	epvector::const_iterator i = seq.begin(), end = seq.end();
	bool coeff_found = false;
	while (i != end) {
		ex t = recombine_pair_to_ex(*i);
		ex c = t.coeff(s, n);
		if (!c.is_zero()) {
			coeffseq.push_back(c);
			coeff_found = 1;
		} else {
			coeffseq.push_back(t);
		}
		++i;
	}
	if (coeff_found) {
		coeffseq.push_back(overall_coeff);
		return (new mul(coeffseq))->setflag(status_flags::dynallocated);
	}
	
	return _ex0;
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x, x1, x2,... stand for a symbolic variables of type ex and c, c1, c2...
 *  stand for such expressions that contain a plain number.
 *  - *(...,x;0) -> 0
 *  - *(+(x1,x2,...);c) -> *(+(*(x1,c),*(x2,c),...))
 *  - *(x;1) -> x
 *  - *(;c) -> c
 *
 *  @param level cut-off in recursive evaluation */
ex mul::eval(int level) const
{
	std::auto_ptr<epvector> evaled_seqp = evalchildren(level);
	if (evaled_seqp.get()) {
		// do more evaluation later
		return (new mul(evaled_seqp, overall_coeff))->
		           setflag(status_flags::dynallocated);
	}
	
#ifdef DO_GINAC_ASSERT
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		GINAC_ASSERT((!is_exactly_a<mul>(i->rest)) ||
		             (!(ex_to<numeric>(i->coeff).is_integer())));
		GINAC_ASSERT(!(i->is_canonical_numeric()));
		if (is_exactly_a<numeric>(recombine_pair_to_ex(*i)))
		    print(print_tree(std::cerr));
		GINAC_ASSERT(!is_exactly_a<numeric>(recombine_pair_to_ex(*i)));
		/* for paranoia */
		expair p = split_ex_to_pair(recombine_pair_to_ex(*i));
		GINAC_ASSERT(p.rest.is_equal(i->rest));
		GINAC_ASSERT(p.coeff.is_equal(i->coeff));
		/* end paranoia */
		++i;
	}
#endif // def DO_GINAC_ASSERT
	
	if (flags & status_flags::evaluated) {
		GINAC_ASSERT(seq.size()>0);
		GINAC_ASSERT(seq.size()>1 || !overall_coeff.is_equal(_ex1));
		return *this;
	}
	
	size_t seq_size = seq.size();
	if (overall_coeff.is_zero()) {
		// *(...,x;0) -> 0
		// unless an element of seq is infinity
		epvector::const_iterator last = seq.end();
		epvector::const_iterator i = seq.begin();
		for (; i != last; ++i) {
			if (i->rest.info(info_flags::infinity) && 
					is_a<numeric>(i->coeff) && 
					ex_to<numeric>(i->coeff).csgn() == 1) {
				throw(std::runtime_error("indeterminate expression: 0*infinity encountered."));
			}
		}
		return _ex0;
	} else if (seq_size==0) {
		// *(;c) -> c
		return overall_coeff;
	} else if (seq_size==1 && seq[0].rest.info(info_flags::infinity)) {
		if (!ex_to<numeric>(overall_coeff).is_real()) {
			throw(std::domain_error("x*Infinity with non real x encountered."));
		} else if (!is_a<numeric>(seq[0].coeff) || 
				!ex_to<numeric>(seq[0].coeff).is_real()) {
				throw(std::domain_error("power::eval(): pow(Infinity, x) for non real x is not defined."));
		} else if (ex_to<numeric>(seq[0].coeff).csgn() == -1) {
			return _ex0;
		} else if (seq[0].rest.is_equal(UnsignedInfinity)) {
			return UnsignedInfinity;
		}
		bool overall_sign=(ex_to<numeric>(overall_coeff).csgn() == -1);
		if (overall_sign) {
			if (seq[0].rest.is_equal(NegInfinity)) {
				if (ex_to<numeric>(seq[0].coeff).is_even()) {
					return NegInfinity;
				} else
					return Infinity;
			} else
				return NegInfinity;
		} else {
			if (seq[0].rest.is_equal(NegInfinity) && 
				ex_to<numeric>(seq[0].coeff).is_even()) {
					return Infinity;
			} else
				return seq[0].rest;
		}
	} else if (seq_size==1 && overall_coeff.is_equal(_ex1) && \
			!ex_to<numeric>(overall_coeff).is_parent_pos_char()) {
		// *(x;1) -> x
		return recombine_pair_to_ex(*(seq.begin()));
	} else if ((seq_size==1) &&
	           is_exactly_a<add>((*seq.begin()).rest) &&
	           ex_to<numeric>((*seq.begin()).coeff).is_equal(*_num1_p)) {
		// *(+(x,y,...);c) -> +(*(x,c),*(y,c),...) (c numeric(), no powers of +())
		const add & addref = ex_to<add>((*seq.begin()).rest);
		std::auto_ptr<epvector> distrseq(new epvector);
		distrseq->reserve(addref.seq.size());
		epvector::const_iterator i = addref.seq.begin(), end = addref.seq.end();
		while (i != end) {
			distrseq->push_back(addref.combine_pair_with_coeff_to_pair(*i, overall_coeff));
			++i;
		}
		return (new add(distrseq,
		                ex_to<numeric>(addref.overall_coeff).
		                mul_dyn(ex_to<numeric>(overall_coeff)))
		       )->setflag(status_flags::dynallocated | status_flags::evaluated);
	} else if ((seq_size >= 2) && (! (flags & status_flags::expanded))) {
		// Strip the content and the unit part from each term. Thus
		// things like (-x+a)*(3*x-3*a) automagically turn into - 3*(x-a)2

		epvector::const_iterator last = seq.end();
		epvector::const_iterator i = seq.begin();
		epvector::const_iterator j = seq.begin();
		std::auto_ptr<epvector> s(new epvector);
		numeric oc = *_num1_p;
		bool something_changed = false;
		// we use pval to store the previously encountered infinity type
		//  0 means we have 1/infinity
		//  1 means we didn't run into infinity
		ex pval = _ex1;
		ex nval = _ex1;
		while (i!=last) {
			if ((i->rest).info(info_flags::infinity)) {
				if (is_a<numeric>(i->coeff) && 
					ex_to<numeric>(i->coeff).is_real()) {
					if (ex_to<numeric>(i->coeff).csgn() ==
							-1){
						nval = _ex0;
					}
					else if (i->rest.is_equal(NegInfinity) 
							&& ex_to<numeric>(
							    i->coeff).is_even())
						nval = Infinity;
					else
						nval = i->rest;
				} else
					throw(std::domain_error("power::eval(): pow(Infinity, x) for non real x is not defined."));
				if (nval.is_zero()) {
					if (!(pval.is_zero() || 
							pval.is_equal(_ex1)))
						throw(std::runtime_error("indeterminate expression: infinity/infinity encountered."));
					else
						pval = _ex0;
				} else if (nval.is_equal(UnsignedInfinity) ||
						pval.is_equal(UnsignedInfinity))
					pval = UnsignedInfinity;
				else if (pval.is_equal(_ex1))
					pval = nval;
				else if (nval.is_equal(pval))
					pval = Infinity;
				else
					pval = NegInfinity;
			}
			if (likely(! (is_a<add>(i->rest) && i->coeff.is_equal(_ex1)))) {
				// power::eval has such a rule, no need to handle powers here
				++i;
				continue;
			}

			// XXX: What is the best way to check if the polynomial is a primitive? 
			numeric c = i->rest.integer_content();
			const numeric lead_coeff =
				ex_to<numeric>(ex_to<add>(i->rest).seq.begin()->coeff).div(c);
			const bool canonicalizable = lead_coeff.is_integer();

			// XXX: The main variable is chosen in a random way, so this code 
			// does NOT transform the term into the canonical form (thus, in some
			// very unlucky event it can even loop forever). Hopefully the main
			// variable will be the same for all terms in *this
			const bool unit_normal = lead_coeff.is_pos_integer();
			if (likely((c == *_num1_p) && ((! canonicalizable) || unit_normal))) {
				++i;
				continue;
			}

			if (! something_changed) {
				s->reserve(seq_size);
				something_changed = true;
			}

			while ((j!=i) && (j!=last)) {
				s->push_back(*j);
				++j;
			}

			if (! unit_normal)
				c = c.mul(*_num_1_p);

			oc = oc.mul(c);

			// divide add by the number in place to save at least 2 .eval() calls
			const add& addref = ex_to<add>(i->rest);
			add* primitive = new add(addref);
			primitive->setflag(status_flags::dynallocated);
			primitive->clearflag(status_flags::hash_calculated);
			primitive->overall_coeff = ex_to<numeric>(primitive->overall_coeff).div_dyn(c);
			for (epvector::iterator ai = primitive->seq.begin();
					ai != primitive->seq.end(); ++ai)
				ai->coeff = ex_to<numeric>(ai->coeff).div_dyn(c);
			
			s->push_back(expair(*primitive, _ex1));

			++i;
			++j;
		}
		if (!pval.is_equal(_ex1)) {
			if (!ex_to<numeric>(overall_coeff).is_real()) {
				throw(std::domain_error("x*Infinity with non real x encountered."));
			} else if (pval.is_zero() || 
					pval.is_equal(UnsignedInfinity))
				return pval;
			else if (ex_to<numeric>(overall_coeff).csgn() == -1)
				return -pval;
			else
				return pval;
		}
		if (something_changed) {
			while (j!=last) {
				s->push_back(*j);
				++j;
			}
			return (new mul(s, ex_to<numeric>(overall_coeff).mul_dyn(oc))
			       )->setflag(status_flags::dynallocated);
		}
	}

	return this->hold();
}

ex mul::evalf(int level, int prec) const
{
	if (level==1)
		return mul(seq,overall_coeff);
	
	if (level==-max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	std::auto_ptr<epvector> s(new epvector);
	s->reserve(seq.size());

	--level;
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		s->push_back(combine_ex_with_coeff_to_pair(
					i->rest.evalf(level, prec), i->coeff));
		++i;
	}
	return mul(s, overall_coeff.evalf(level, prec));
}

void mul::find_real_imag(ex & rp, ex & ip) const
{
	rp = overall_coeff.real_part();
	ip = overall_coeff.imag_part();
	for (epvector::const_iterator i=seq.begin(); i!=seq.end(); ++i) {
		ex factor = recombine_pair_to_ex(*i);
		ex new_rp = factor.real_part();
		ex new_ip = factor.imag_part();
		if(new_ip.is_zero()) {
			rp *= new_rp;
			ip *= new_rp;
		} else {
			ex temp = rp*new_rp - ip*new_ip;
			ip = ip*new_rp + rp*new_ip;
			rp = temp;
		}
	}
	rp = rp.expand();
	ip = ip.expand();
}

ex mul::real_part() const
{
	ex rp, ip;
	find_real_imag(rp, ip);
	return rp;
}

ex mul::imag_part() const
{
	ex rp, ip;
	find_real_imag(rp, ip);
	return ip;
}

ex mul::evalm() const
{
	// numeric*matrix
	if (seq.size() == 1 && seq[0].coeff.is_equal(_ex1)
	 && is_a<matrix>(seq[0].rest))
		return ex_to<matrix>(seq[0].rest).mul(ex_to<numeric>(overall_coeff));

	// Evaluate children first, look whether there are any matrices at all
	// (there can be either no matrices or one matrix; if there were more
	// than one matrix, it would be a non-commutative product)
	std::auto_ptr<epvector> s(new epvector);
	s->reserve(seq.size());

	bool have_matrix = false;
	epvector::iterator the_matrix;

	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		const ex &m = recombine_pair_to_ex(*i).evalm();
		s->push_back(split_ex_to_pair(m));
		if (is_a<matrix>(m)) {
			have_matrix = true;
			the_matrix = s->end() - 1;
		}
		++i;
	}

	if (have_matrix) {

		// The product contained a matrix. We will multiply all other factors
		// into that matrix.
		matrix m = ex_to<matrix>(the_matrix->rest);
		s->erase(the_matrix);
		ex scalar = (new mul(s, overall_coeff))->setflag(status_flags::dynallocated);
		return m.mul_scalar(scalar);

	} else
		return (new mul(s, overall_coeff))->setflag(status_flags::dynallocated);
}

ex mul::eval_ncmul(const exvector & v) const
{
	if (seq.empty())
		return inherited::eval_ncmul(v);

	// Find first noncommutative element and call its eval_ncmul()
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		if (i->rest.return_type() == return_types::noncommutative)
			return i->rest.eval_ncmul(v);
		++i;
	}
	return inherited::eval_ncmul(v);
}

bool tryfactsubs(const ex & origfactor, const ex & patternfactor, int & nummatches, lst & repls)
{	
	ex origbase;
	int origexponent;
	int origexpsign;

	if (is_exactly_a<power>(origfactor) && origfactor.op(1).info(info_flags::integer)) {
		origbase = origfactor.op(0);
		int expon = ex_to<numeric>(origfactor.op(1)).to_int();
		origexponent = expon > 0 ? expon : -expon;
		origexpsign = expon > 0 ? 1 : -1;
	} else {
		origbase = origfactor;
		origexponent = 1;
		origexpsign = 1;
	}

	ex patternbase;
	int patternexponent;
	int patternexpsign;

	if (is_exactly_a<power>(patternfactor) && patternfactor.op(1).info(info_flags::integer)) {
		patternbase = patternfactor.op(0);
		int expon = ex_to<numeric>(patternfactor.op(1)).to_int();
		patternexponent = expon > 0 ? expon : -expon;
		patternexpsign = expon > 0 ? 1 : -1;
	} else {
		patternbase = patternfactor;
		patternexponent = 1;
		patternexpsign = 1;
	}

	lst saverepls = repls;
	if (origexponent < patternexponent || origexpsign != patternexpsign || !origbase.match(patternbase,saverepls))
		return false;
	repls = saverepls;

	int newnummatches = origexponent / patternexponent;
	if (newnummatches < nummatches)
		nummatches = newnummatches;
	return true;
}

/** Checks wheter e matches to the pattern pat and the (possibly to be updated)
  * list of replacements repls. This matching is in the sense of algebraic
  * substitutions. Matching starts with pat.op(factor) of the pattern because
  * the factors before this one have already been matched. The (possibly
  * updated) number of matches is in nummatches. subsed[i] is true for factors
  * that already have been replaced by previous substitutions and matched[i]
  * is true for factors that have been matched by the current match.
  */
bool algebraic_match_mul_with_mul(const mul &e, const ex &pat, lst &repls,
		int factor, int &nummatches, const std::vector<bool> &subsed,
		std::vector<bool> &matched)
{
	if (factor == pat.nops())
		return true;

	for (size_t i=0; i<e.nops(); ++i) {
		if(subsed[i] || matched[i])
			continue;
		lst newrepls = repls;
		int newnummatches = nummatches;
		if (tryfactsubs(e.op(i), pat.op(factor), newnummatches, newrepls)) {
			matched[i] = true;
			if (algebraic_match_mul_with_mul(e, pat, newrepls, factor+1,
					newnummatches, subsed, matched)) {
				repls = newrepls;
				nummatches = newnummatches;
				return true;
			}
			else
				matched[i] = false;
		}
	}

	return false;
}

bool mul::has(const ex & pattern, unsigned options) const
{
	if(!(options&has_options::algebraic))
		return basic::has(pattern,options);
	if(is_a<mul>(pattern)) {
		lst repls;
		int nummatches = std::numeric_limits<int>::max();
		std::vector<bool> subsed(seq.size(), false);
		std::vector<bool> matched(seq.size(), false);
		if(algebraic_match_mul_with_mul(*this, pattern, repls, 0, nummatches,
				subsed, matched))
			return true;
	}
	return basic::has(pattern, options);
}

ex mul::algebraic_subs_mul(const exmap & m, unsigned options) const
{	
	std::vector<bool> subsed(seq.size(), false);
	exvector subsresult(seq.size());
	ex divide_by = 1;
	ex multiply_by = 1;

	for (exmap::const_iterator it = m.begin(); it != m.end(); ++it) {

		if (is_exactly_a<mul>(it->first)) {
retry1:
			int nummatches = std::numeric_limits<int>::max();
			std::vector<bool> currsubsed(seq.size(), false);
			lst repls;
			
			if(!algebraic_match_mul_with_mul(*this, it->first, repls, 0, nummatches, subsed, currsubsed))
				continue;

			for (size_t j=0; j<subsed.size(); j++)
				if (currsubsed[j])
					subsed[j] = true;
			ex subsed_pattern
				= it->first.subs(ex(repls), subs_options::no_pattern);
			divide_by *= power(subsed_pattern, nummatches);
			ex subsed_result
				= it->second.subs(ex(repls), subs_options::no_pattern);
			multiply_by *= power(subsed_result, nummatches);
			goto retry1;

		} else {

			for (size_t j=0; j<this->nops(); j++) {
				int nummatches = std::numeric_limits<int>::max();
				lst repls;
				if (!subsed[j] && tryfactsubs(op(j), it->first, nummatches, repls)){
					subsed[j] = true;
					ex subsed_pattern
						= it->first.subs(ex(repls), subs_options::no_pattern);
					divide_by *= power(subsed_pattern, nummatches);
					ex subsed_result
						= it->second.subs(ex(repls), subs_options::no_pattern);
					multiply_by *= power(subsed_result, nummatches);
				}
			}
		}
	}

	bool subsfound = false;
	for (size_t i=0; i<subsed.size(); i++) {
		if (subsed[i]) {
			subsfound = true;
			break;
		}
	}
	if (!subsfound)
		return subs_one_level(m, options | subs_options::algebraic);

	return ((*this)/divide_by)*multiply_by;
}

// protected

/** Implementation of ex::diff() for a product.  It applies the product rule.
 *  @see ex::diff */
ex mul::derivative(const symbol & s) const
{
	size_t num = seq.size();
	exvector addseq;
	addseq.reserve(num);
	
	// D(a*b*c) = D(a)*b*c + a*D(b)*c + a*b*D(c)
	epvector mulseq = seq;
	epvector::const_iterator i = seq.begin(), end = seq.end();
	epvector::iterator i2 = mulseq.begin();
	while (i != end) {
		expair ep = split_ex_to_pair(power(i->rest, i->coeff - _ex1) *
		                             i->rest.diff(s));
		ep.swap(*i2);
		addseq.push_back((new mul(mulseq, overall_coeff * i->coeff))->setflag(status_flags::dynallocated));
		ep.swap(*i2);
		++i; ++i2;
	}
	return (new add(addseq))->setflag(status_flags::dynallocated);
}

double mul::total_degree() const
{
	if (flags & status_flags::tdegree_calculated) {
		return tdegree;
	}
	numeric tdeg = calc_total_degree();
	if (tdeg.is_real())
		tdegree = tdeg.to_double();
	else
		tdegree = std::sqrt(std::pow(tdeg.real().to_double(), 2) + 
				std::pow(tdeg.imag().to_double(), 2));
	setflag(status_flags::tdegree_calculated);
	return tdegree;
}

int mul::compare(const basic& other) const
{
	static const tinfo_t function_id = find_tinfo_key("function");
	static const tinfo_t fderivative_id = find_tinfo_key("fderivative");
	static const tinfo_t pow_id = find_tinfo_key("power");
	static const tinfo_t symbol_id = find_tinfo_key("symbol");
	const tinfo_t typeid_this = tinfo();
	const tinfo_t typeid_other = other.tinfo();
	if (typeid_this==typeid_other) {
		GINAC_ASSERT(typeid(*this)==typeid(other));
		return compare_same_type(other);
	} else if (typeid_other == pow_id) {
		return compare_pow(static_cast<const power&>(other));
	} else if (typeid_other == symbol_id) {
		return compare_symbol(static_cast<const symbol&>(other));
	} else if (typeid_other == function_id ||
			typeid_other == fderivative_id) {
		return -1;
	} else {
		return (typeid_this<typeid_other ? -1 : 1);
	}
}

int mul::compare_symbol(const symbol &other) const
{
	int cmpval;
	double tdeg;
	tdeg = total_degree();
	if (tdeg == 1) {
		cmpval = seq[0].rest.compare(other);
		if (cmpval != 0) {
			return cmpval;
		}
		cmpval = _ex1.compare(seq[0].coeff);
		if (cmpval != 0) {
			return cmpval;
		}
		return -1;
	}
	return tdeg > 1 ? -1 : 1;
}

// compare this to a pow object
// first we compare degrees
// if equal we compare the first item in the sequence to the base in other
int mul::compare_pow(const power &other) const
{
	double my_deg = total_degree();
	double odeg;
	numeric oexp;
	int cmpval = 0;
	if (is_a<numeric>(other.exponent)) {
		numeric oexp = ex_to<numeric>(other.exponent);
		if (oexp.is_real()) {
			odeg = oexp.to_double();
		} else {
			odeg = std::sqrt(std::pow(oexp.real().to_double(), 2) + 
					std::pow(oexp.imag().to_double(), 2));
		}
		if (odeg != my_deg)
			return my_deg < odeg ? 1 : -1;
	} else {
		cmpval = seq[0].coeff.compare(other.exponent);
		if (cmpval != 0)
			return cmpval;
	}
	cmpval = seq[0].rest.compare(other.basis);
	if (cmpval != 0) {
		return cmpval;
	}
	if (seq.size() == 1 && overall_coeff.is_equal(_ex_1))
		return 0;
	return 1;
}


int mul::compare_same_type(const basic & other) const
{
	int cmpval;

	const mul &o = static_cast<const mul &>(other);

	// compare total degrees
	double deg1 = total_degree();
	double deg2 = o.total_degree();
	if (deg1 != deg2)
		return deg1 < deg2 ? 1 : -1;

	// compare each item in this product to correnponding element in other
	epvector::const_iterator cit1 = seq.begin();
	epvector::const_iterator cit2 = o.seq.begin();
	epvector::const_iterator last1 = seq.end();
	epvector::const_iterator last2 = o.seq.end();

	for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
		cmpval = (*cit1).compare(*cit2);
		if (cmpval != 0)
			return cmpval;
	}

	// compare sizes
	if (cit1 != last1) 
		return 1;
	else if (cit2 != last2)
		return -1;

	return 0;
}

unsigned mul::return_type() const
{
	if (seq.empty()) {
		// mul without factors: should not happen, but commutates
		return return_types::commutative;
	}
	
	bool all_commutative = true;
	epvector::const_iterator noncommutative_element; // point to first found nc element
	
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		unsigned rt = i->rest.return_type();
		if (rt == return_types::noncommutative_composite)
			return rt; // one ncc -> mul also ncc
		if ((rt == return_types::noncommutative) && (all_commutative)) {
			// first nc element found, remember position
			noncommutative_element = i;
			all_commutative = false;
		}
		if ((rt == return_types::noncommutative) && (!all_commutative)) {
			// another nc element found, compare type_infos
			if (noncommutative_element->rest.return_type_tinfo() != i->rest.return_type_tinfo()) {
					// different types -> mul is ncc
					return return_types::noncommutative_composite;
			}
		}
		++i;
	}
	// all factors checked
	return all_commutative ? return_types::commutative : return_types::noncommutative;
}
   
tinfo_t mul::return_type_tinfo() const
{
	if (seq.empty())
		return this;  // mul without factors: should not happen
	
	// return type_info of first noncommutative element
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		if (i->rest.return_type() == return_types::noncommutative)
			return i->rest.return_type_tinfo();
		++i;
	}
	// no noncommutative element found, should not happen
	return this;
}

ex mul::thisexpairseq(const epvector & v, const ex & oc, bool do_index_renaming) const
{
	return (new mul(v, oc, do_index_renaming))->setflag(status_flags::dynallocated);
}

ex mul::thisexpairseq(std::auto_ptr<epvector> vp, const ex & oc, bool do_index_renaming) const
{
	return (new mul(vp, oc, do_index_renaming))->setflag(status_flags::dynallocated);
}

expair mul::split_ex_to_pair(const ex & e) const
{
	if (is_exactly_a<power>(e)) {
		const power & powerref = ex_to<power>(e);
		if (is_exactly_a<numeric>(powerref.exponent))
			return expair(powerref.basis,powerref.exponent);
	}
	return expair(e,_ex1);
}
	
expair mul::combine_ex_with_coeff_to_pair(const ex & e,
                                          const ex & c) const
{
	// to avoid duplication of power simplification rules,
	// we create a temporary power object
	// otherwise it would be hard to correctly evaluate
	// expression like (4^(1/3))^(3/2)
	if (c.is_equal(_ex1))
		return split_ex_to_pair(e);

	return split_ex_to_pair(power(e,c));
}
	
expair mul::combine_pair_with_coeff_to_pair(const expair & p,
                                            const ex & c) const
{
	// to avoid duplication of power simplification rules,
	// we create a temporary power object
	// otherwise it would be hard to correctly evaluate
	// expression like (4^(1/3))^(3/2)
	if (c.is_equal(_ex1))
		return p;

	return split_ex_to_pair(power(recombine_pair_to_ex(p),c));
}
	
ex mul::recombine_pair_to_ex(const expair & p) const
{
	if (ex_to<numeric>(p.coeff).is_equal(*_num1_p)) 
		return p.rest;
	else
		return (new power(p.rest,p.coeff))->setflag(status_flags::dynallocated);
}

bool mul::expair_needs_further_processing(epp it)
{
	if (is_exactly_a<mul>(it->rest) &&
		ex_to<numeric>(it->coeff).is_integer()) {
		// combined pair is product with integer power -> expand it
		*it = split_ex_to_pair(recombine_pair_to_ex(*it));
		return true;
	}
	if (is_exactly_a<numeric>(it->rest)) {
		expair ep = split_ex_to_pair(recombine_pair_to_ex(*it));
		if (!ep.is_equal(*it)) {
			// combined pair is a numeric power which can be simplified
			*it = ep;
			return true;
		}
		if (it->coeff.is_equal(_ex1)) {
			// combined pair has coeff 1 and must be moved to the end
			return true;
		}
	}
	return false;
}       

ex mul::default_overall_coeff() const
{
	return _ex1;
}

void mul::combine_overall_coeff(const ex & c)
{
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	overall_coeff = ex_to<numeric>(overall_coeff).mul_dyn(ex_to<numeric>(c));
}

void mul::combine_overall_coeff(const ex & c1, const ex & c2)
{
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c1));
	GINAC_ASSERT(is_exactly_a<numeric>(c2));
	overall_coeff = ex_to<numeric>(overall_coeff).mul_dyn(ex_to<numeric>(c1).power(ex_to<numeric>(c2)));
}

bool mul::can_make_flat(const expair & p) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));
	// this assertion will probably fail somewhere
	// it would require a more careful make_flat, obeying the power laws
	// probably should return true only if p.coeff is integer
	return ex_to<numeric>(p.coeff).is_equal(*_num1_p);
}

bool mul::can_be_further_expanded(const ex & e)
{
	if (is_exactly_a<mul>(e)) {
		for (epvector::const_iterator cit = ex_to<mul>(e).seq.begin(); cit != ex_to<mul>(e).seq.end(); ++cit) {
			if (is_exactly_a<add>(cit->rest) && cit->coeff.info(info_flags::posint))
				return true;
		}
	} else if (is_exactly_a<power>(e)) {
		if (is_exactly_a<add>(e.op(0)) && e.op(1).info(info_flags::posint))
			return true;
	}
	return false;
}

ex mul::expand(unsigned options) const
{
	{
	// trivial case: expanding the monomial (~ 30% of all calls)
		epvector::const_iterator i = seq.begin(), seq_end = seq.end();
		while ((i != seq.end()) &&  is_a<symbol>(i->rest) && i->coeff.info(info_flags::integer))
			++i;
		if (i == seq_end) {
			setflag(status_flags::expanded);
			return *this;
		}
	}

	// do not rename indices if the object has no indices at all
	if ((!(options & expand_options::expand_rename_idx)) && 
			this->info(info_flags::has_indices))
		options |= expand_options::expand_rename_idx;

	const bool skip_idx_rename = !(options & expand_options::expand_rename_idx);

	// First, expand the children
	std::auto_ptr<epvector> expanded_seqp = expandchildren(options);
	const epvector & expanded_seq = (expanded_seqp.get() ? *expanded_seqp : seq);

	// Now, look for all the factors that are sums and multiply each one out
	// with the next one that is found while collecting the factors which are
	// not sums
	ex last_expanded = _ex1;

	epvector non_adds;
	non_adds.reserve(expanded_seq.size());

	for (epvector::const_iterator cit = expanded_seq.begin(); cit != expanded_seq.end(); ++cit) {
		if (is_exactly_a<add>(cit->rest) &&
			(cit->coeff.is_equal(_ex1))) {
			if (is_exactly_a<add>(last_expanded)) {

				// Expand a product of two sums, aggressive version.
				// Caring for the overall coefficients in separate loops can
				// sometimes give a performance gain of up to 15%!

				const int sizedifference = ex_to<add>(last_expanded).seq.size()-ex_to<add>(cit->rest).seq.size();
				// add2 is for the inner loop and should be the bigger of the two sums
				// in the presence of asymptotically good sorting:
				const add& add1 = (sizedifference<0 ? ex_to<add>(last_expanded) : ex_to<add>(cit->rest));
				const add& add2 = (sizedifference<0 ? ex_to<add>(cit->rest) : ex_to<add>(last_expanded));
				const epvector::const_iterator add1begin = add1.seq.begin();
				const epvector::const_iterator add1end   = add1.seq.end();
				const epvector::const_iterator add2begin = add2.seq.begin();
				const epvector::const_iterator add2end   = add2.seq.end();
				epvector distrseq;
				distrseq.reserve(add1.seq.size()+add2.seq.size());

				// Multiply add2 with the overall coefficient of add1 and append it to distrseq:
				if (!add1.overall_coeff.is_zero()) {
					if (add1.overall_coeff.is_equal(_ex1))
						distrseq.insert(distrseq.end(),add2begin,add2end);
					else
						for (epvector::const_iterator i=add2begin; i!=add2end; ++i)
							distrseq.push_back(expair(i->rest, ex_to<numeric>(i->coeff).mul_dyn(ex_to<numeric>(add1.overall_coeff))));
				}

				// Multiply add1 with the overall coefficient of add2 and append it to distrseq:
				if (!add2.overall_coeff.is_zero()) {
					if (add2.overall_coeff.is_equal(_ex1))
						distrseq.insert(distrseq.end(),add1begin,add1end);
					else
						for (epvector::const_iterator i=add1begin; i!=add1end; ++i)
							distrseq.push_back(expair(i->rest, ex_to<numeric>(i->coeff).mul_dyn(ex_to<numeric>(add2.overall_coeff))));
				}

				// Compute the new overall coefficient and put it together:
				ex tmp_accu = (new add(distrseq, add1.overall_coeff*add2.overall_coeff))->setflag(status_flags::dynallocated);

				exvector add1_dummy_indices, add2_dummy_indices, add_indices;
				lst dummy_subs;

				if (!skip_idx_rename) {
					for (epvector::const_iterator i=add1begin; i!=add1end; ++i) {
						add_indices = get_all_dummy_indices_safely(i->rest);
						add1_dummy_indices.insert(add1_dummy_indices.end(), add_indices.begin(), add_indices.end());
					}
					for (epvector::const_iterator i=add2begin; i!=add2end; ++i) {
						add_indices = get_all_dummy_indices_safely(i->rest);
						add2_dummy_indices.insert(add2_dummy_indices.end(), add_indices.begin(), add_indices.end());
					}

					sort(add1_dummy_indices.begin(), add1_dummy_indices.end(), ex_is_less());
					sort(add2_dummy_indices.begin(), add2_dummy_indices.end(), ex_is_less());
					dummy_subs = rename_dummy_indices_uniquely(add1_dummy_indices, add2_dummy_indices);
				}

				// Multiply explicitly all non-numeric terms of add1 and add2:
				for (epvector::const_iterator i2=add2begin; i2!=add2end; ++i2) {
					// We really have to combine terms here in order to compactify
					// the result.  Otherwise it would become waayy tooo bigg.
					numeric oc(*_num0_p);
					epvector distrseq2;
					distrseq2.reserve(add1.seq.size());
					const ex i2_new = (skip_idx_rename || (dummy_subs.op(0).nops() == 0) ?
							i2->rest :
							i2->rest.subs(ex_to<lst>(dummy_subs.op(0)), 
								ex_to<lst>(dummy_subs.op(1)), subs_options::no_pattern));
					for (epvector::const_iterator i1=add1begin; i1!=add1end; ++i1) {
						// Don't push_back expairs which might have a rest that evaluates to a numeric,
						// since that would violate an invariant of expairseq:
						const ex rest = (new mul(i1->rest, i2_new))->setflag(status_flags::dynallocated);
						if (is_exactly_a<numeric>(rest)) {
							oc += ex_to<numeric>(rest).mul(ex_to<numeric>(i1->coeff).mul(ex_to<numeric>(i2->coeff)));
						} else {
							distrseq2.push_back(expair(rest, ex_to<numeric>(i1->coeff).mul_dyn(ex_to<numeric>(i2->coeff))));
						}
					}
					tmp_accu += (new add(distrseq2, oc))->setflag(status_flags::dynallocated);
				} 
				last_expanded = tmp_accu;
			} else {
				if (!last_expanded.is_equal(_ex1))
					non_adds.push_back(split_ex_to_pair(last_expanded));
				last_expanded = cit->rest;
			}

		} else {
			non_adds.push_back(*cit);
		}
	}

	// Now the only remaining thing to do is to multiply the factors which
	// were not sums into the "last_expanded" sum
	if (is_exactly_a<add>(last_expanded)) {
		size_t n = last_expanded.nops();
		exvector distrseq;
		distrseq.reserve(n);
		exvector va;
		if (! skip_idx_rename) {
			va = get_all_dummy_indices_safely(mul(non_adds));
			sort(va.begin(), va.end(), ex_is_less());
		}

		for (size_t i=0; i<n; ++i) {
			epvector factors = non_adds;
			if (skip_idx_rename)
				factors.push_back(split_ex_to_pair(last_expanded.op(i)));
			else
				factors.push_back(split_ex_to_pair(rename_dummy_indices_uniquely(va, last_expanded.op(i))));
			ex term = (new mul(factors, overall_coeff))->setflag(status_flags::dynallocated);
			if (can_be_further_expanded(term)) {
				distrseq.push_back(term.expand());
			} else {
				if (options == 0)
					ex_to<basic>(term).setflag(status_flags::expanded);
				distrseq.push_back(term);
			}
		}

		return ((new add(distrseq))->
		        setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0)));
	}

	non_adds.push_back(split_ex_to_pair(last_expanded));
	ex result = (new mul(non_adds, overall_coeff))->setflag(status_flags::dynallocated);
	if (can_be_further_expanded(result)) {
		return result.expand();
	} else {
		if (options == 0)
			ex_to<basic>(result).setflag(status_flags::expanded);
		return result;
	}
}

  
//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////


/** Member-wise expand the expairs representing this sequence.  This must be
 *  overridden from expairseq::expandchildren() and done iteratively in order
 *  to allow for early cancallations and thus safe memory.
 *
 *  @see mul::expand()
 *  @return pointer to epvector containing expanded representation or zero
 *  pointer, if sequence is unchanged. */
std::auto_ptr<epvector> mul::expandchildren(unsigned options) const
{
	const epvector::const_iterator last = seq.end();
	epvector::const_iterator cit = seq.begin();
	while (cit!=last) {
		const ex & factor = recombine_pair_to_ex(*cit);
		const ex & expanded_factor = factor.expand(options);
		if (!are_ex_trivially_equal(factor,expanded_factor)) {
			
			// something changed, copy seq, eval and return it
			std::auto_ptr<epvector> s(new epvector);
			s->reserve(seq.size());
			
			// copy parts of seq which are known not to have changed
			epvector::const_iterator cit2 = seq.begin();
			while (cit2!=cit) {
				s->push_back(*cit2);
				++cit2;
			}

			// copy first changed element
			s->push_back(split_ex_to_pair(expanded_factor));
			++cit2;

			// copy rest
			while (cit2!=last) {
				s->push_back(split_ex_to_pair(recombine_pair_to_ex(*cit2).expand(options)));
				++cit2;
			}
			return s;
		}
		++cit;
	}
	
	return std::auto_ptr<epvector>(0); // nothing has changed
}

} // namespace GiNaC
