/** @file ncmul.cpp
 *
 *  Implementation of GiNaC's non-commutative products of expressions. */

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

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "ncmul.h"
#include "ex.h"
#include "add.h"
#include "mul.h"
#include "clifford.h"
#include "matrix.h"
#include "archive.h"
#include "indexed.h"
#include "utils.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(ncmul, exprseq,
  print_func<print_context>(&ncmul::do_print).
  print_func<print_tree>(&ncmul::do_print_tree).
  print_func<print_csrc>(&ncmul::do_print_csrc).
  print_func<print_python_repr>(&ncmul::do_print_csrc))


//////////
// default constructor
//////////

ncmul::ncmul()
{
	tinfo_key = &ncmul::tinfo_static;
}

//////////
// other constructors
//////////

// public

ncmul::ncmul(const ex & lh, const ex & rh) : inherited(lh,rh)
{
	tinfo_key = &ncmul::tinfo_static;
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3) : inherited(f1,f2,f3)
{
	tinfo_key = &ncmul::tinfo_static;
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3,
             const ex & f4) : inherited(f1,f2,f3,f4)
{
	tinfo_key = &ncmul::tinfo_static;
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3,
             const ex & f4, const ex & f5) : inherited(f1,f2,f3,f4,f5)
{
	tinfo_key = &ncmul::tinfo_static;
}

ncmul::ncmul(const ex & f1, const ex & f2, const ex & f3,
             const ex & f4, const ex & f5, const ex & f6) : inherited(f1,f2,f3,f4,f5,f6)
{
	tinfo_key = &ncmul::tinfo_static;
}

ncmul::ncmul(const exvector & v, bool discardable) : inherited(v,discardable)
{
	tinfo_key = &ncmul::tinfo_static;
}

ncmul::ncmul(std::auto_ptr<exvector> vp) : inherited(vp)
{
	tinfo_key = &ncmul::tinfo_static;
}

//////////
// archiving
//////////

DEFAULT_ARCHIVING(ncmul)
	
//////////
// functions overriding virtual functions from base classes
//////////

// public

void ncmul::do_print(const print_context & c, unsigned level) const
{
	printseq(c, "(", '*', ")", precedence(), level);
}

void ncmul::do_print_csrc(const print_context & c, unsigned level) const
{
	c.s << class_name();
	printseq(c, "(", ',', ")", precedence(), precedence());
}

bool ncmul::info(unsigned inf) const
{
	return inherited::info(inf);
}

typedef std::vector<int> intvector;

ex ncmul::expand(unsigned options) const
{
	// First, expand the children
	std::auto_ptr<exvector> vp = expandchildren(options);
	const exvector &expanded_seq = vp.get() ? *vp : this->seq;
	
	// Now, look for all the factors that are sums and remember their
	// position and number of terms.
	intvector positions_of_adds(expanded_seq.size());
	intvector number_of_add_operands(expanded_seq.size());

	size_t number_of_adds = 0;
	size_t number_of_expanded_terms = 1;

	size_t current_position = 0;
	exvector::const_iterator last = expanded_seq.end();
	for (exvector::const_iterator cit=expanded_seq.begin(); cit!=last; ++cit) {
		if (is_exactly_a<add>(*cit)) {
			positions_of_adds[number_of_adds] = current_position;
			size_t num_ops = cit->nops();
			number_of_add_operands[number_of_adds] = num_ops;
			number_of_expanded_terms *= num_ops;
			number_of_adds++;
		}
		++current_position;
	}

	// If there are no sums, we are done
	if (number_of_adds == 0) {
		if (vp.get())
			return (new ncmul(vp))->
			        setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
		else
			return *this;
	}

	// Now, form all possible products of the terms of the sums with the
	// remaining factors, and add them together
	exvector distrseq;
	distrseq.reserve(number_of_expanded_terms);

	intvector k(number_of_adds);

	/* Rename indices in the static members of the product */
	exvector expanded_seq_mod;
	size_t j = 0;
	exvector va;

	for (size_t i=0; i<expanded_seq.size(); i++) {
		if (i == positions_of_adds[j]) {
			expanded_seq_mod.push_back(_ex1);
			j++;
		} else {
			expanded_seq_mod.push_back(rename_dummy_indices_uniquely(va, expanded_seq[i], true));
		}
	}

	while (true) {
		exvector term = expanded_seq_mod;
		for (size_t i=0; i<number_of_adds; i++) {
			term[positions_of_adds[i]] = rename_dummy_indices_uniquely(va, expanded_seq[positions_of_adds[i]].op(k[i]), true);
		}

		distrseq.push_back((new ncmul(term, true))->
		                    setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0)));

		// increment k[]
		int l = number_of_adds-1;
		while ((l>=0) && ((++k[l]) >= number_of_add_operands[l])) {
			k[l] = 0;
			l--;
		}
		if (l<0)
			break;
	}

	return (new add(distrseq))->
	        setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
}

int ncmul::degree(const ex & s) const
{
	if (is_equal(ex_to<basic>(s)))
		return 1;

	// Sum up degrees of factors
	int deg_sum = 0;
	exvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		deg_sum += i->degree(s);
		++i;
	}
	return deg_sum;
}

int ncmul::ldegree(const ex & s) const
{
	if (is_equal(ex_to<basic>(s)))
		return 1;

	// Sum up degrees of factors
	int deg_sum = 0;
	exvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		deg_sum += i->degree(s);
		++i;
	}
	return deg_sum;
}

ex ncmul::coeff(const ex & s, int n) const
{
	if (is_equal(ex_to<basic>(s)))
		return n==1 ? _ex1 : _ex0;

	exvector coeffseq;
	coeffseq.reserve(seq.size());

	if (n == 0) {
		// product of individual coeffs
		// if a non-zero power of s is found, the resulting product will be 0
		exvector::const_iterator it=seq.begin();
		while (it!=seq.end()) {
			coeffseq.push_back((*it).coeff(s,n));
			++it;
		}
		return (new ncmul(coeffseq,1))->setflag(status_flags::dynallocated);
	}
		 
	exvector::const_iterator i = seq.begin(), end = seq.end();
	bool coeff_found = false;
	while (i != end) {
		ex c = i->coeff(s,n);
		if (c.is_zero()) {
			coeffseq.push_back(*i);
		} else {
			coeffseq.push_back(c);
			coeff_found = true;
		}
		++i;
	}

	if (coeff_found) return (new ncmul(coeffseq,1))->setflag(status_flags::dynallocated);
	
	return _ex0;
}

size_t ncmul::count_factors(const ex & e) const
{
	if ((is_exactly_a<mul>(e)&&(e.return_type()!=return_types::commutative))||
		(is_exactly_a<ncmul>(e))) {
		size_t factors=0;
		for (size_t i=0; i<e.nops(); i++)
			factors += count_factors(e.op(i));
		
		return factors;
	}
	return 1;
}
		
void ncmul::append_factors(exvector & v, const ex & e) const
{
	if ((is_exactly_a<mul>(e)&&(e.return_type()!=return_types::commutative))||
		(is_exactly_a<ncmul>(e))) {
		for (size_t i=0; i<e.nops(); i++)
			append_factors(v, e.op(i));
	} else 
		v.push_back(e);
}

typedef std::vector<unsigned> unsignedvector;
typedef std::vector<exvector> exvectorvector;

/** Perform automatic term rewriting rules in this class.  In the following
 *  x, x1, x2,... stand for a symbolic variables of type ex and c, c1, c2...
 *  stand for such expressions that contain a plain number.
 *  - ncmul(...,*(x1,x2),...,ncmul(x3,x4),...) -> ncmul(...,x1,x2,...,x3,x4,...)  (associativity)
 *  - ncmul(x) -> x
 *  - ncmul() -> 1
 *  - ncmul(...,c1,...,c2,...) -> *(c1,c2,ncmul(...))  (pull out commutative elements)
 *  - ncmul(x1,y1,x2,y2) -> *(ncmul(x1,x2),ncmul(y1,y2))  (collect elements of same type)
 *  - ncmul(x1,x2,x3,...) -> x::eval_ncmul(x1,x2,x3,...)
 *
 *  @param level cut-off in recursive evaluation */
ex ncmul::eval(int level) const
{
	// The following additional rule would be nice, but produces a recursion,
	// which must be trapped by introducing a flag that the sub-ncmuls()
	// are already evaluated (maybe later...)
	//                  ncmul(x1,x2,...,X,y1,y2,...) ->
	//                      ncmul(ncmul(x1,x2,...),X,ncmul(y1,y2,...)
	//                      (X noncommutative_composite)

	if ((level==1) && (flags & status_flags::evaluated)) {
		return *this;
	}

	exvector evaledseq=evalchildren(level);

	// ncmul(...,*(x1,x2),...,ncmul(x3,x4),...) ->
	//     ncmul(...,x1,x2,...,x3,x4,...)  (associativity)
	size_t factors = 0;
	exvector::const_iterator cit = evaledseq.begin(), citend = evaledseq.end();
	while (cit != citend)
		factors += count_factors(*cit++);
	
	exvector assocseq;
	assocseq.reserve(factors);
	cit = evaledseq.begin();
	make_flat_inserter mf(evaledseq, true);
	while (cit != citend)
	{	ex factor = mf.handle_factor(*(cit++), 1);
		append_factors(assocseq, factor);
	}
	
	// ncmul(x) -> x
	if (assocseq.size()==1) return *(seq.begin());
 
	// ncmul() -> 1
	if (assocseq.empty()) return _ex1;

	// determine return types
	unsignedvector rettypes;
	rettypes.reserve(assocseq.size());
	size_t i = 0;
	size_t count_commutative=0;
	size_t count_noncommutative=0;
	size_t count_noncommutative_composite=0;
	cit = assocseq.begin(); citend = assocseq.end();
	while (cit != citend) {
		switch (rettypes[i] = cit->return_type()) {
		case return_types::commutative:
			count_commutative++;
			break;
		case return_types::noncommutative:
			count_noncommutative++;
			break;
		case return_types::noncommutative_composite:
			count_noncommutative_composite++;
			break;
		default:
			throw(std::logic_error("ncmul::eval(): invalid return type"));
		}
		++i; ++cit;
	}
	GINAC_ASSERT(count_commutative+count_noncommutative+count_noncommutative_composite==assocseq.size());

	// ncmul(...,c1,...,c2,...) ->
	//     *(c1,c2,ncmul(...)) (pull out commutative elements)
	if (count_commutative!=0) {
		exvector commutativeseq;
		commutativeseq.reserve(count_commutative+1);
		exvector noncommutativeseq;
		noncommutativeseq.reserve(assocseq.size()-count_commutative);
		size_t num = assocseq.size();
		for (size_t i=0; i<num; ++i) {
			if (rettypes[i]==return_types::commutative)
				commutativeseq.push_back(assocseq[i]);
			else
				noncommutativeseq.push_back(assocseq[i]);
		}
		commutativeseq.push_back((new ncmul(noncommutativeseq,1))->setflag(status_flags::dynallocated));
		return (new mul(commutativeseq))->setflag(status_flags::dynallocated);
	}
		
	// ncmul(x1,y1,x2,y2) -> *(ncmul(x1,x2),ncmul(y1,y2))
	//     (collect elements of same type)

	if (count_noncommutative_composite==0) {
		// there are neither commutative nor noncommutative_composite
		// elements in assocseq
		GINAC_ASSERT(count_commutative==0);

		size_t assoc_num = assocseq.size();
		exvectorvector evv;
		std::vector<tinfo_t> rttinfos;
		evv.reserve(assoc_num);
		rttinfos.reserve(assoc_num);

		cit = assocseq.begin(), citend = assocseq.end();
		while (cit != citend) {
			tinfo_t ti = cit->return_type_tinfo();
			size_t rtt_num = rttinfos.size();
			// search type in vector of known types
			for (i=0; i<rtt_num; ++i) {
				if(ti == rttinfos[i]) {
					evv[i].push_back(*cit);
					break;
				}
			}
			if (i >= rtt_num) {
				// new type
				rttinfos.push_back(ti);
				evv.push_back(exvector());
				(evv.end()-1)->reserve(assoc_num);
				(evv.end()-1)->push_back(*cit);
			}
			++cit;
		}

		size_t evv_num = evv.size();
#ifdef DO_GINAC_ASSERT
		GINAC_ASSERT(evv_num == rttinfos.size());
		GINAC_ASSERT(evv_num > 0);
		size_t s=0;
		for (i=0; i<evv_num; ++i)
			s += evv[i].size();
		GINAC_ASSERT(s == assoc_num);
#endif // def DO_GINAC_ASSERT
		
		// if all elements are of same type, simplify the string
		if (evv_num == 1) {
			return evv[0][0].eval_ncmul(evv[0]);
		}
		
		exvector splitseq;
		splitseq.reserve(evv_num);
		for (i=0; i<evv_num; ++i)
			splitseq.push_back((new ncmul(evv[i]))->setflag(status_flags::dynallocated));
		
		return (new mul(splitseq))->setflag(status_flags::dynallocated);
	}
	
	return (new ncmul(assocseq))->setflag(status_flags::dynallocated |
										  status_flags::evaluated);
}

ex ncmul::evalm() const
{
	// Evaluate children first
	std::auto_ptr<exvector> s(new exvector);
	s->reserve(seq.size());
	exvector::const_iterator it = seq.begin(), itend = seq.end();
	while (it != itend) {
		s->push_back(it->evalm());
		it++;
	}

	// If there are only matrices, simply multiply them
	it = s->begin(); itend = s->end();
	if (is_a<matrix>(*it)) {
		matrix prod(ex_to<matrix>(*it));
		it++;
		while (it != itend) {
			if (!is_a<matrix>(*it))
				goto no_matrix;
			prod = prod.mul(ex_to<matrix>(*it));
			it++;
		}
		return prod;
	}

no_matrix:
	return (new ncmul(s))->setflag(status_flags::dynallocated);
}

ex ncmul::thiscontainer(const exvector & v) const
{
	return (new ncmul(v))->setflag(status_flags::dynallocated);
}

ex ncmul::thiscontainer(std::auto_ptr<exvector> vp) const
{
	return (new ncmul(vp))->setflag(status_flags::dynallocated);
}

ex ncmul::conjugate() const
{
	if (return_type() != return_types::noncommutative) {
		return exprseq::conjugate();
	}

	if (!is_clifford_tinfo(return_type_tinfo())) {
		return exprseq::conjugate();
	}

	exvector ev;
	ev.reserve(nops());
	for (const_iterator i=end(); i!=begin();) {
		--i;
		ev.push_back(i->conjugate());
	}
	return (new ncmul(ev, true))->setflag(status_flags::dynallocated).eval();
}

ex ncmul::real_part() const
{
	return basic::real_part();
}

ex ncmul::imag_part() const
{
	return basic::imag_part();
}

// protected

/** Implementation of ex::diff() for a non-commutative product. It applies
 *  the product rule.
 *  @see ex::diff */
ex ncmul::derivative(const symbol & s) const
{
	size_t num = seq.size();
	exvector addseq;
	addseq.reserve(num);
	
	// D(a*b*c) = D(a)*b*c + a*D(b)*c + a*b*D(c)
	exvector ncmulseq = seq;
	for (size_t i=0; i<num; ++i) {
		ex e = seq[i].diff(s);
		e.swap(ncmulseq[i]);
		addseq.push_back((new ncmul(ncmulseq))->setflag(status_flags::dynallocated));
		e.swap(ncmulseq[i]);
	}
	return (new add(addseq))->setflag(status_flags::dynallocated);
}

int ncmul::compare_same_type(const basic & other) const
{
	return inherited::compare_same_type(other);
}

unsigned ncmul::return_type() const
{
	if (seq.empty())
		return return_types::commutative;

	bool all_commutative = true;
	exvector::const_iterator noncommutative_element; // point to first found nc element

	exvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		unsigned rt = i->return_type();
		if (rt == return_types::noncommutative_composite)
			return rt; // one ncc -> mul also ncc
		if ((rt == return_types::noncommutative) && (all_commutative)) {
			// first nc element found, remember position
			noncommutative_element = i;
			all_commutative = false;
		}
		if ((rt == return_types::noncommutative) && (!all_commutative)) {
			// another nc element found, compare type_infos
			if(noncommutative_element->return_type_tinfo() != i->return_type_tinfo())
					return return_types::noncommutative_composite;
		}
		++i;
	}
	// all factors checked
	GINAC_ASSERT(!all_commutative); // not all factors should commutate, because this is a ncmul();
	return all_commutative ? return_types::commutative : return_types::noncommutative;
}
   
tinfo_t ncmul::return_type_tinfo() const
{
	if (seq.empty())
		return this;

	// return type_info of first noncommutative element
	exvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		if (i->return_type() == return_types::noncommutative)
			return i->return_type_tinfo();
		++i;
	}

	// no noncommutative element found, should not happen
	return this;
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

std::auto_ptr<exvector> ncmul::expandchildren(unsigned options) const
{
	const_iterator cit = this->seq.begin(), end = this->seq.end();
	while (cit != end) {
		const ex & expanded_ex = cit->expand(options);
		if (!are_ex_trivially_equal(*cit, expanded_ex)) {

			// copy first part of seq which hasn't changed
			std::auto_ptr<exvector> s(new exvector(this->seq.begin(), cit));
			reserve(*s, this->seq.size());

			// insert changed element
			s->push_back(expanded_ex);
			++cit;

			// copy rest
			while (cit != end) {
				s->push_back(cit->expand(options));
				++cit;
			}

			return s;
		}

		++cit;
	}

	return std::auto_ptr<exvector>(0); // nothing has changed
}

const exvector & ncmul::get_factors() const
{
	return seq;
}

//////////
// friend functions
//////////

ex reeval_ncmul(const exvector & v)
{
	return (new ncmul(v))->setflag(status_flags::dynallocated);
}

ex hold_ncmul(const exvector & v)
{
	if (v.empty())
		return _ex1;
	else if (v.size() == 1)
		return v[0];
	else
		return (new ncmul(v))->setflag(status_flags::dynallocated |
		                               status_flags::evaluated);
}

} // namespace GiNaC
