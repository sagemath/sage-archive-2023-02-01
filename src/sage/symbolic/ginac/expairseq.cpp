/** @file expairseq.cpp
 *
 *  Implementation of sequences of expression pairs. */

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
#include <algorithm>
#include <string>
#include <stdexcept>
#include <iterator>

#include "expairseq.h"
#include "lst.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "relational.h"
#include "wildcard.h"
#include "archive.h"
#include "operators.h"
#include "utils.h"
#include "indexed.h"

#if EXPAIRSEQ_USE_HASHTAB
#include <cmath>
#endif // EXPAIRSEQ_USE_HASHTAB

namespace GiNaC {

	
GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(expairseq, basic,
  print_func<print_context>(&expairseq::do_print).
  print_func<print_tree>(&expairseq::do_print_tree))


//////////
// helper classes
//////////

class epp_is_less
{
public:
	bool operator()(const epp &lh, const epp &rh) const
	{
		return (*lh).is_less(*rh);
	}
};

//////////
// default constructor
//////////

// public

expairseq::expairseq() : inherited(&expairseq::tinfo_static)
#if EXPAIRSEQ_USE_HASHTAB
                                                   , hashtabsize(0)
#endif // EXPAIRSEQ_USE_HASHTAB
{}

// protected

#if 0
/** For use by copy ctor and assignment operator. */
void expairseq::copy(const expairseq &other)
{
	seq = other.seq;
	overall_coeff = other.overall_coeff;
#if EXPAIRSEQ_USE_HASHTAB
	// copy hashtab
	hashtabsize = other.hashtabsize;
	if (hashtabsize!=0) {
		hashmask = other.hashmask;
		hashtab.resize(hashtabsize);
		epvector::const_iterator osb = other.seq.begin();
		for (unsigned i=0; i<hashtabsize; ++i) {
			hashtab[i].clear();
			for (epplist::const_iterator cit=other.hashtab[i].begin();
			     cit!=other.hashtab[i].end(); ++cit) {
				hashtab[i].push_back(seq.begin()+((*cit)-osb));
			}
		}
	} else {
		hashtab.clear();
	}
#endif // EXPAIRSEQ_USE_HASHTAB
}
#endif

//////////
// other constructors
//////////

expairseq::expairseq(const ex &lh, const ex &rh) : inherited(&expairseq::tinfo_static)
{
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(const exvector &v) : inherited(&expairseq::tinfo_static)
{
	construct_from_exvector(v);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(const epvector &v, const ex &oc, bool do_index_renaming)
  : inherited(&expairseq::tinfo_static), overall_coeff(oc)
{
	GINAC_ASSERT(is_a<numeric>(oc));
	construct_from_epvector(v, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(std::auto_ptr<epvector> vp, const ex &oc, bool do_index_renaming)
  : inherited(&expairseq::tinfo_static), overall_coeff(oc)
{
	GINAC_ASSERT(vp.get()!=0);
	GINAC_ASSERT(is_a<numeric>(oc));
	construct_from_epvector(*vp, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

//////////
// archiving
//////////

expairseq::expairseq(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
#if EXPAIRSEQ_USE_HASHTAB
	, hashtabsize(0)
#endif
{
	archive_node::archive_node_cit first = n.find_first("rest");
	archive_node::archive_node_cit last = n.find_last("coeff");
	++last;
	seq.reserve((last-first)/2);

	for (archive_node::archive_node_cit loc = first; loc < last;) {
		ex rest;
		ex coeff;
		n.find_ex_by_loc(loc++, rest, sym_lst);
		n.find_ex_by_loc(loc++, coeff, sym_lst);
		seq.push_back(expair(rest, coeff));
	}

	n.find_ex("overall_coeff", overall_coeff, sym_lst);

	canonicalize();
	GINAC_ASSERT(is_canonical());
}

void expairseq::archive(archive_node &n) const
{
	inherited::archive(n);
	epvector::const_iterator i = seq.begin(), iend = seq.end();
	while (i != iend) {
		n.add_ex("rest", i->rest);
		n.add_ex("coeff", i->coeff);
		++i;
	}
	n.add_ex("overall_coeff", overall_coeff);
}

DEFAULT_UNARCHIVE(expairseq)

//////////
// functions overriding virtual functions from base classes
//////////

// public

void expairseq::do_print(const print_context & c, unsigned level) const
{
	c.s << "[[";
	printseq(c, ',', precedence(), level);
	c.s << "]]";
}

void expairseq::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", nops=" << nops()
	    << std::endl;
	size_t num = seq.size();
	for (size_t i=0; i<num; ++i) {
		seq[i].rest.print(c, level + c.delta_indent);
		seq[i].coeff.print(c, level + c.delta_indent);
		if (i != num - 1)
			c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl;
	}
	if (!overall_coeff.is_equal(default_overall_coeff())) {
		c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl
		    << std::string(level + c.delta_indent, ' ') << "overall_coeff" << std::endl;
		overall_coeff.print(c, level + c.delta_indent);
	}
	c.s << std::string(level + c.delta_indent,' ') << "=====" << std::endl;
#if EXPAIRSEQ_USE_HASHTAB
	c.s << std::string(level + c.delta_indent,' ')
	    << "hashtab size " << hashtabsize << std::endl;
	if (hashtabsize == 0) return;
#define MAXCOUNT 5
	unsigned count[MAXCOUNT+1];
	for (int i=0; i<MAXCOUNT+1; ++i)
		count[i] = 0;
	unsigned this_bin_fill;
	unsigned cum_fill_sq = 0;
	unsigned cum_fill = 0;
	for (unsigned i=0; i<hashtabsize; ++i) {
		this_bin_fill = 0;
		if (hashtab[i].size() > 0) {
			c.s << std::string(level + c.delta_indent, ' ')
			    << "bin " << i << " with entries ";
			for (epplist::const_iterator it=hashtab[i].begin();
			     it!=hashtab[i].end(); ++it) {
				c.s << *it-seq.begin() << " ";
				++this_bin_fill;
			}
			c.s << std::endl;
			cum_fill += this_bin_fill;
			cum_fill_sq += this_bin_fill*this_bin_fill;
		}
		if (this_bin_fill<MAXCOUNT)
			++count[this_bin_fill];
		else
			++count[MAXCOUNT];
	}
	unsigned fact = 1;
	double cum_prob = 0;
	double lambda = (1.0*seq.size()) / hashtabsize;
	for (int k=0; k<MAXCOUNT; ++k) {
		if (k>0)
			fact *= k;
		double prob = std::pow(lambda,k)/fact * std::exp(-lambda);
		cum_prob += prob;
		c.s << std::string(level + c.delta_indent, ' ') << "bins with " << k << " entries: "
		    << int(1000.0*count[k]/hashtabsize)/10.0 << "% (expected: "
		    << int(prob*1000)/10.0 << ")" << std::endl;
	}
	c.s << std::string(level + c.delta_indent, ' ') << "bins with more entries: "
	    << int(1000.0*count[MAXCOUNT]/hashtabsize)/10.0 << "% (expected: "
	    << int((1-cum_prob)*1000)/10.0 << ")" << std::endl;

	c.s << std::string(level + c.delta_indent, ' ') << "variance: "
	    << 1.0/hashtabsize*cum_fill_sq-(1.0/hashtabsize*cum_fill)*(1.0/hashtabsize*cum_fill)
	    << std::endl;
	c.s << std::string(level + c.delta_indent, ' ') << "average fill: "
	    << (1.0*cum_fill)/hashtabsize
	    << " (should be equal to " << (1.0*seq.size())/hashtabsize << ")" << std::endl;
#endif // EXPAIRSEQ_USE_HASHTAB
}

bool expairseq::info(unsigned inf) const
{
	switch(inf) {
		case info_flags::expanded:
			return (flags & status_flags::expanded);
		case info_flags::has_indices: {
			if (flags & status_flags::has_indices)
				return true;
			else if (flags & status_flags::has_no_indices)
				return false;
			for (epvector::const_iterator i = seq.begin(); i != seq.end(); ++i) {
				if (i->rest.info(info_flags::has_indices)) {
					this->setflag(status_flags::has_indices);
					this->clearflag(status_flags::has_no_indices);
					return true;
				}
			}
			this->clearflag(status_flags::has_indices);
			this->setflag(status_flags::has_no_indices);
			return false;
		}
	}
	return inherited::info(inf);
}

size_t expairseq::nops() const
{
	if (overall_coeff.is_equal(default_overall_coeff()))
		return seq.size();
	else
		return seq.size()+1;
}

ex expairseq::op(size_t i) const
{
	if (i < seq.size())
		return recombine_pair_to_ex(seq[i]);
	GINAC_ASSERT(!overall_coeff.is_equal(default_overall_coeff()));
	return overall_coeff;
}

ex expairseq::map(map_function &f) const
{
	std::auto_ptr<epvector> v(new epvector);
	v->reserve(seq.size()+1);

	epvector::const_iterator cit = seq.begin(), last = seq.end();
	while (cit != last) {
		v->push_back(split_ex_to_pair(f(recombine_pair_to_ex(*cit))));
		++cit;
	}

	if (overall_coeff.is_equal(default_overall_coeff()))
		return thisexpairseq(v, default_overall_coeff(), true);
	else {
		ex newcoeff = f(overall_coeff);
		if(is_a<numeric>(newcoeff))
			return thisexpairseq(v, newcoeff, true);
		else {
			v->push_back(split_ex_to_pair(newcoeff));
			return thisexpairseq(v, default_overall_coeff(), true);
		}
	}
}

/** Perform coefficient-wise automatic term rewriting rules in this class. */
ex expairseq::eval(int level) const
{
	if ((level==1) && (flags &status_flags::evaluated))
		return *this;
	
	std::auto_ptr<epvector> vp = evalchildren(level);
	if (vp.get() == 0)
		return this->hold();
	
	return (new expairseq(vp, overall_coeff))->setflag(status_flags::dynallocated | status_flags::evaluated);
}

epvector* conjugateepvector(const epvector&epv)
{
	epvector *newepv = 0;
	for (epvector::const_iterator i=epv.begin(); i!=epv.end(); ++i) {
		if(newepv) {
			newepv->push_back(i->conjugate());
			continue;
		}
		expair x = i->conjugate();
		if (x.is_equal(*i)) {
			continue;
		}
		newepv = new epvector;
		newepv->reserve(epv.size());
		for (epvector::const_iterator j=epv.begin(); j!=i; ++j) {
			newepv->push_back(*j);
		}
		newepv->push_back(x);
	}
	return newepv;
}

ex expairseq::conjugate() const
{
	epvector* newepv = conjugateepvector(seq);
	ex x = overall_coeff.conjugate();
	if (!newepv && are_ex_trivially_equal(x, overall_coeff)) {
		return *this;
	}
	ex result = thisexpairseq(newepv ? *newepv : seq, x);
	if (newepv) {
		delete newepv;
	}
	return result;
}

bool expairseq::is_polynomial(const ex & var) const
{
	if (!is_exactly_a<add>(*this) && !is_exactly_a<mul>(*this))
		return basic::is_polynomial(var);
	for (epvector::const_iterator i=seq.begin(); i!=seq.end(); ++i) {
		if (!(i->rest).is_polynomial(var))
			return false;
	}
	return true;
}

bool expairseq::match(const ex & pattern, lst & repl_lst) const
{
	// This differs from basic::match() because we want "a+b+c+d" to
	// match "d+*+b" with "*" being "a+c", and we want to honor commutativity

	if (this->tinfo() == ex_to<basic>(pattern).tinfo()) {

		// Check whether global wildcard (one that matches the "rest of the
		// expression", like "*" above) is present
		bool has_global_wildcard = false;
		ex global_wildcard;
		for (size_t i=0; i<pattern.nops(); i++) {
			if (is_exactly_a<wildcard>(pattern.op(i))) {
				has_global_wildcard = true;
				global_wildcard = pattern.op(i);
				break;
			}
		}

		// Unfortunately, this is an O(N^2) operation because we can't
		// sort the pattern in a useful way...

		// Chop into terms
		exvector ops;
		ops.reserve(nops());
		for (size_t i=0; i<nops(); i++)
			ops.push_back(op(i));

		// Now, for every term of the pattern, look for a matching term in
		// the expression and remove the match
		for (size_t i=0; i<pattern.nops(); i++) {
			ex p = pattern.op(i);
			if (has_global_wildcard && p.is_equal(global_wildcard))
				continue;
			exvector::iterator it = ops.begin(), itend = ops.end();
			while (it != itend) {
				lst::const_iterator last_el = repl_lst.end();
				--last_el;
				if (it->match(p, repl_lst)) {
					ops.erase(it);
					goto found;
				}
				while(true) {
					lst::const_iterator next_el = last_el;
					++next_el;
					if(next_el == repl_lst.end())
						break;
					else
						repl_lst.remove_last();
				}
				++it;
			}
			return false; // no match found
found:		;
		}

		if (has_global_wildcard) {

			// Assign all the remaining terms to the global wildcard (unless
			// it has already been matched before, in which case the matches
			// must be equal)
			size_t num = ops.size();
			std::auto_ptr<epvector> vp(new epvector);
			vp->reserve(num);
			for (size_t i=0; i<num; i++)
				vp->push_back(split_ex_to_pair(ops[i]));
			ex rest = thisexpairseq(vp, default_overall_coeff());
			for (lst::const_iterator it = repl_lst.begin(); it != repl_lst.end(); ++it) {
				if (it->op(0).is_equal(global_wildcard))
					return rest.is_equal(it->op(1));
			}
			repl_lst.append(global_wildcard == rest);
			return true;

		} else {

			// No global wildcard, then the match fails if there are any
			// unmatched terms left
			return ops.empty();
		}
	}
	return inherited::match(pattern, repl_lst);
}

ex expairseq::subs(const exmap & m, unsigned options) const
{
	std::auto_ptr<epvector> vp = subschildren(m, options);
	if (vp.get())
		return ex_to<basic>(thisexpairseq(vp, overall_coeff, true));
	else if ((options & subs_options::algebraic) && is_exactly_a<mul>(*this))
		return static_cast<const mul *>(this)->algebraic_subs_mul(m, options);
	else
		return subs_one_level(m, options);
}

// protected

int expairseq::compare_same_type(const basic &other) const
{
	GINAC_ASSERT(is_a<expairseq>(other));
	const expairseq &o = static_cast<const expairseq &>(other);
	
	int cmpval;
	
	// compare number of elements
	if (seq.size() != o.seq.size())
		return (seq.size()<o.seq.size()) ? -1 : 1;
	
	// compare overall_coeff
	cmpval = overall_coeff.compare(o.overall_coeff);
	if (cmpval!=0)
		return cmpval;
	
#if EXPAIRSEQ_USE_HASHTAB
	GINAC_ASSERT(hashtabsize==o.hashtabsize);
	if (hashtabsize==0) {
#endif // EXPAIRSEQ_USE_HASHTAB
		epvector::const_iterator cit1 = seq.begin();
		epvector::const_iterator cit2 = o.seq.begin();
		epvector::const_iterator last1 = seq.end();
		epvector::const_iterator last2 = o.seq.end();
		
		for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
			cmpval = (*cit1).compare(*cit2);
			if (cmpval!=0) return cmpval;
		}
		
		GINAC_ASSERT(cit1==last1);
		GINAC_ASSERT(cit2==last2);
		
		return 0;
#if EXPAIRSEQ_USE_HASHTAB
	}
	
	// compare number of elements in each hashtab entry
	for (unsigned i=0; i<hashtabsize; ++i) {
		unsigned cursize=hashtab[i].size();
		if (cursize != o.hashtab[i].size())
			return (cursize < o.hashtab[i].size()) ? -1 : 1;
	}
	
	// compare individual (sorted) hashtab entries
	for (unsigned i=0; i<hashtabsize; ++i) {
		unsigned sz = hashtab[i].size();
		if (sz>0) {
			const epplist &eppl1 = hashtab[i];
			const epplist &eppl2 = o.hashtab[i];
			epplist::const_iterator it1 = eppl1.begin();
			epplist::const_iterator it2 = eppl2.begin();
			while (it1!=eppl1.end()) {
				cmpval = (*(*it1)).compare(*(*it2));
				if (cmpval!=0)
					return cmpval;
				++it1;
				++it2;
			}
		}
	}
	
	return 0; // equal
#endif // EXPAIRSEQ_USE_HASHTAB
}

bool expairseq::is_equal_same_type(const basic &other) const
{
	const expairseq &o = static_cast<const expairseq &>(other);
	
	// compare number of elements
	if (seq.size()!=o.seq.size())
		return false;
	
	// compare overall_coeff
	if (!overall_coeff.is_equal(o.overall_coeff))
		return false;
	
#if EXPAIRSEQ_USE_HASHTAB
	// compare number of elements in each hashtab entry
	if (hashtabsize!=o.hashtabsize) {
		std::cout << "this:" << std::endl;
		print(print_tree(std::cout));
		std::cout << "other:" << std::endl;
		other.print(print_tree(std::cout));
	}
		
	GINAC_ASSERT(hashtabsize==o.hashtabsize);
	
	if (hashtabsize==0) {
#endif // EXPAIRSEQ_USE_HASHTAB
		epvector::const_iterator cit1 = seq.begin();
		epvector::const_iterator cit2 = o.seq.begin();
		epvector::const_iterator last1 = seq.end();
		
		while (cit1!=last1) {
			if (!(*cit1).is_equal(*cit2)) return false;
			++cit1;
			++cit2;
		}
		
		return true;
#if EXPAIRSEQ_USE_HASHTAB
	}
	
	for (unsigned i=0; i<hashtabsize; ++i) {
		if (hashtab[i].size() != o.hashtab[i].size())
			return false;
	}

	// compare individual sorted hashtab entries
	for (unsigned i=0; i<hashtabsize; ++i) {
		unsigned sz = hashtab[i].size();
		if (sz>0) {
			const epplist &eppl1 = hashtab[i];
			const epplist &eppl2 = o.hashtab[i];
			epplist::const_iterator it1 = eppl1.begin();
			epplist::const_iterator it2 = eppl2.begin();
			while (it1!=eppl1.end()) {
				if (!(*(*it1)).is_equal(*(*it2))) return false;
				++it1;
				++it2;
			}
		}
	}
	
	return true;
#endif // EXPAIRSEQ_USE_HASHTAB
}

unsigned expairseq::return_type() const
{
	return return_types::noncommutative_composite;
}

unsigned expairseq::calchash() const
{
	unsigned v = golden_ratio_hash((p_int)this->tinfo());
	epvector::const_iterator i = seq.begin();
	const epvector::const_iterator end = seq.end();
	while (i != end) {
		v ^= i->rest.gethash();
#if !EXPAIRSEQ_USE_HASHTAB
		// rotation spoils commutativity!
		v = rotate_left(v);
		v ^= i->coeff.gethash();
#endif // !EXPAIRSEQ_USE_HASHTAB
		++i;
	}

	v ^= overall_coeff.gethash();

	// store calculated hash value only if object is already evaluated
	if (flags &status_flags::evaluated) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}
	
	return v;
}

ex expairseq::expand(unsigned options) const
{
	std::auto_ptr<epvector> vp = expandchildren(options);
	if (vp.get())
		return thisexpairseq(vp, overall_coeff);
	else {
		// The terms have not changed, so it is safe to declare this expanded
		return (options == 0) ? setflag(status_flags::expanded) : *this;
	}
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// protected

/** Create an object of this type.
 *  This method works similar to a constructor.  It is useful because expairseq
 *  has (at least) two possible different semantics but we want to inherit
 *  methods thus avoiding code duplication.  Sometimes a method in expairseq
 *  has to create a new one of the same semantics, which cannot be done by a
 *  ctor because the name (add, mul,...) is unknown on the expaiseq level.  In
 *  order for this trick to work a derived class must of course override this
 *  definition. */
ex expairseq::thisexpairseq(const epvector &v, const ex &oc, bool do_index_renaming) const
{
	return expairseq(v, oc, do_index_renaming);
}

ex expairseq::thisexpairseq(std::auto_ptr<epvector> vp, const ex &oc, bool do_index_renaming) const
{
	return expairseq(vp, oc, do_index_renaming);
}

void expairseq::printpair(const print_context & c, const expair & p, unsigned upper_precedence) const
{
	c.s << "[[";
	p.rest.print(c, precedence());
	c.s << ",";
	p.coeff.print(c, precedence());
	c.s << "]]";
}

void expairseq::printseq(const print_context & c, char delim,
                         unsigned this_precedence,
                         unsigned upper_precedence) const
{
	if (this_precedence <= upper_precedence)
		c.s << "(";
	epvector::const_iterator it, it_last = seq.end() - 1;
	for (it=seq.begin(); it!=it_last; ++it) {
		printpair(c, *it, this_precedence);
		c.s << delim;
	}
	printpair(c, *it, this_precedence);
	if (!overall_coeff.is_equal(default_overall_coeff())) {
		c.s << delim;
		overall_coeff.print(c, this_precedence);
	}
	
	if (this_precedence <= upper_precedence)
		c.s << ")";
}


/** Form an expair from an ex, using the corresponding semantics.
 *  @see expairseq::recombine_pair_to_ex() */
expair expairseq::split_ex_to_pair(const ex &e) const
{
	return expair(e,_ex1);
}


expair expairseq::combine_ex_with_coeff_to_pair(const ex &e,
                                                const ex &c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	
	return expair(e,c);
}


expair expairseq::combine_pair_with_coeff_to_pair(const expair &p,
                                                  const ex &c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	
	return expair(p.rest,ex_to<numeric>(p.coeff).mul_dyn(ex_to<numeric>(c)));
}


/** Form an ex out of an expair, using the corresponding semantics.
 *  @see expairseq::split_ex_to_pair() */
ex expairseq::recombine_pair_to_ex(const expair &p) const
{
	return lst(p.rest,p.coeff);
}

bool expairseq::expair_needs_further_processing(epp it)
{
#if EXPAIRSEQ_USE_HASHTAB
	//#  error "FIXME: expair_needs_further_processing not yet implemented for hashtabs, sorry. A.F."
#endif // EXPAIRSEQ_USE_HASHTAB
	return false;
}

ex expairseq::default_overall_coeff() const
{
	return _ex0;
}

void expairseq::combine_overall_coeff(const ex &c)
{
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	overall_coeff = ex_to<numeric>(overall_coeff).add_dyn(ex_to<numeric>(c));
}

void expairseq::combine_overall_coeff(const ex &c1, const ex &c2)
{
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c1));
	GINAC_ASSERT(is_exactly_a<numeric>(c2));
	overall_coeff = ex_to<numeric>(overall_coeff).
	                add_dyn(ex_to<numeric>(c1).mul(ex_to<numeric>(c2)));
}

bool expairseq::can_make_flat(const expair &p) const
{
	return true;
}


//////////
// non-virtual functions in this class
//////////

void expairseq::construct_from_2_ex_via_exvector(const ex &lh, const ex &rh)
{
	exvector v;
	v.reserve(2);
	v.push_back(lh);
	v.push_back(rh);
	construct_from_exvector(v);
#if EXPAIRSEQ_USE_HASHTAB
	GINAC_ASSERT((hashtabsize==0)||(hashtabsize>=minhashtabsize));
	GINAC_ASSERT(hashtabsize==calc_hashtabsize(seq.size()));
#endif // EXPAIRSEQ_USE_HASHTAB
}

void expairseq::construct_from_2_ex(const ex &lh, const ex &rh)
{
	if (ex_to<basic>(lh).tinfo()==this->tinfo()) {
		if (ex_to<basic>(rh).tinfo()==this->tinfo()) {
#if EXPAIRSEQ_USE_HASHTAB
			unsigned totalsize = ex_to<expairseq>(lh).seq.size() +
			                     ex_to<expairseq>(rh).seq.size();
			if (calc_hashtabsize(totalsize)!=0) {
				construct_from_2_ex_via_exvector(lh,rh);
			} else {
#endif // EXPAIRSEQ_USE_HASHTAB
				if (is_a<mul>(lh) && lh.info(info_flags::has_indices) && 
					rh.info(info_flags::has_indices)) {
					ex newrh=rename_dummy_indices_uniquely(lh, rh);
					construct_from_2_expairseq(ex_to<expairseq>(lh),
					                           ex_to<expairseq>(newrh));
				}
				else
					construct_from_2_expairseq(ex_to<expairseq>(lh),
					                           ex_to<expairseq>(rh));
#if EXPAIRSEQ_USE_HASHTAB
			}
#endif // EXPAIRSEQ_USE_HASHTAB
			return;
		} else {
#if EXPAIRSEQ_USE_HASHTAB
			unsigned totalsize = ex_to<expairseq>(lh).seq.size()+1;
			if (calc_hashtabsize(totalsize)!=0) {
				construct_from_2_ex_via_exvector(lh, rh);
			} else {
#endif // EXPAIRSEQ_USE_HASHTAB
				construct_from_expairseq_ex(ex_to<expairseq>(lh), rh);
#if EXPAIRSEQ_USE_HASHTAB
			}
#endif // EXPAIRSEQ_USE_HASHTAB
			return;
		}
	} else if (ex_to<basic>(rh).tinfo()==this->tinfo()) {
#if EXPAIRSEQ_USE_HASHTAB
		unsigned totalsize=ex_to<expairseq>(rh).seq.size()+1;
		if (calc_hashtabsize(totalsize)!=0) {
			construct_from_2_ex_via_exvector(lh,rh);
		} else {
#endif // EXPAIRSEQ_USE_HASHTAB
			construct_from_expairseq_ex(ex_to<expairseq>(rh),lh);
#if EXPAIRSEQ_USE_HASHTAB
		}
#endif // EXPAIRSEQ_USE_HASHTAB
		return;
	}
	
#if EXPAIRSEQ_USE_HASHTAB
	if (calc_hashtabsize(2)!=0) {
		construct_from_2_ex_via_exvector(lh,rh);
		return;
	}
	hashtabsize = 0;
#endif // EXPAIRSEQ_USE_HASHTAB
	
	if (is_exactly_a<numeric>(lh)) {
		if (is_exactly_a<numeric>(rh)) {
			combine_overall_coeff(lh);
			combine_overall_coeff(rh);
		} else {
			combine_overall_coeff(lh);
			seq.push_back(split_ex_to_pair(rh));
		}
	} else {
		if (is_exactly_a<numeric>(rh)) {
			combine_overall_coeff(rh);
			seq.push_back(split_ex_to_pair(lh));
		} else {
			expair p1 = split_ex_to_pair(lh);
			expair p2 = split_ex_to_pair(rh);
			
			int cmpval = p1.rest.compare(p2.rest);
			if (cmpval==0) {
				p1.coeff = ex_to<numeric>(p1.coeff).add_dyn(ex_to<numeric>(p2.coeff));
				if (!ex_to<numeric>(p1.coeff).is_zero()) {
					// no further processing is necessary, since this
					// one element will usually be recombined in eval()
					seq.push_back(p1);
				}
			} else {
				seq.reserve(2);
				if (cmpval<0) {
					seq.push_back(p1);
					seq.push_back(p2);
				} else {
					seq.push_back(p2);
					seq.push_back(p1);
				}
			}
		}
	}
}

void expairseq::construct_from_2_expairseq(const expairseq &s1,
										   const expairseq &s2)
{
	combine_overall_coeff(s1.overall_coeff);
	combine_overall_coeff(s2.overall_coeff);

	epvector::const_iterator first1 = s1.seq.begin();
	epvector::const_iterator last1 = s1.seq.end();
	epvector::const_iterator first2 = s2.seq.begin();
	epvector::const_iterator last2 = s2.seq.end();

	seq.reserve(s1.seq.size()+s2.seq.size());

	bool needs_further_processing=false;
	
	while (first1!=last1 && first2!=last2) {
		int cmpval = (*first1).rest.compare((*first2).rest);

		if (cmpval==0) {
			// combine terms
			const numeric &newcoeff = ex_to<numeric>(first1->coeff).
			                           add(ex_to<numeric>(first2->coeff));
			if (!newcoeff.is_zero()) {
				seq.push_back(expair(first1->rest,newcoeff));
				if (expair_needs_further_processing(seq.end()-1)) {
					needs_further_processing = true;
				}
			}
			++first1;
			++first2;
		} else if (cmpval<0) {
			seq.push_back(*first1);
			++first1;
		} else {
			seq.push_back(*first2);
			++first2;
		}
	}
	
	while (first1!=last1) {
		seq.push_back(*first1);
		++first1;
	}
	while (first2!=last2) {
		seq.push_back(*first2);
		++first2;
	}
	
	if (needs_further_processing) {
		epvector v = seq;
		seq.clear();
		construct_from_epvector(v);
	}
}

void expairseq::construct_from_expairseq_ex(const expairseq &s,
											const ex &e)
{
	combine_overall_coeff(s.overall_coeff);
	if (is_exactly_a<numeric>(e)) {
		combine_overall_coeff(e);
		seq = s.seq;
		return;
	}
	
	epvector::const_iterator first = s.seq.begin();
	epvector::const_iterator last = s.seq.end();
	expair p = split_ex_to_pair(e);
	
	seq.reserve(s.seq.size()+1);
	bool p_pushed = false;
	
	bool needs_further_processing=false;
	
	// merge p into s.seq
	while (first!=last) {
		int cmpval = (*first).rest.compare(p.rest);
		if (cmpval==0) {
			// combine terms
			const numeric &newcoeff = ex_to<numeric>(first->coeff).
			                           add(ex_to<numeric>(p.coeff));
			if (!newcoeff.is_zero()) {
				seq.push_back(expair(first->rest,newcoeff));
				if (expair_needs_further_processing(seq.end()-1))
					needs_further_processing = true;
			}
			++first;
			p_pushed = true;
			break;
		} else if (cmpval<0) {
			seq.push_back(*first);
			++first;
		} else {
			seq.push_back(p);
			p_pushed = true;
			break;
		}
	}
	
	if (p_pushed) {
		// while loop exited because p was pushed, now push rest of s.seq
		while (first!=last) {
			seq.push_back(*first);
			++first;
		}
	} else {
		// while loop exited because s.seq was pushed, now push p
		seq.push_back(p);
	}

	if (needs_further_processing) {
		epvector v = seq;
		seq.clear();
		construct_from_epvector(v);
	}
}

void expairseq::construct_from_exvector(const exvector &v)
{
	// simplifications: +(a,+(b,c),d) -> +(a,b,c,d) (associativity)
	//                  +(d,b,c,a) -> +(a,b,c,d) (canonicalization)
	//                  +(...,x,*(x,c1),*(x,c2)) -> +(...,*(x,1+c1+c2)) (c1, c2 numeric())
	//                  (same for (+,*) -> (*,^)

	make_flat(v);
#if EXPAIRSEQ_USE_HASHTAB
	combine_same_terms();
#else
	canonicalize();
	combine_same_terms_sorted_seq();
#endif // EXPAIRSEQ_USE_HASHTAB
}

void expairseq::construct_from_epvector(const epvector &v, bool do_index_renaming)
{
	// simplifications: +(a,+(b,c),d) -> +(a,b,c,d) (associativity)
	//                  +(d,b,c,a) -> +(a,b,c,d) (canonicalization)
	//                  +(...,x,*(x,c1),*(x,c2)) -> +(...,*(x,1+c1+c2)) (c1, c2 numeric())
	//                  (same for (+,*) -> (*,^)

	make_flat(v, do_index_renaming);
#if EXPAIRSEQ_USE_HASHTAB
	combine_same_terms();
#else
	canonicalize();
	combine_same_terms_sorted_seq();
#endif // EXPAIRSEQ_USE_HASHTAB
}

/** Combine this expairseq with argument exvector.
 *  It cares for associativity as well as for special handling of numerics. */
void expairseq::make_flat(const exvector &v)
{
	exvector::const_iterator cit;
	
	// count number of operands which are of same expairseq derived type
	// and their cumulative number of operands
	int nexpairseqs = 0;
	int noperands = 0;
	bool do_idx_rename = false;
	
	cit = v.begin();
	while (cit!=v.end()) {
		if (ex_to<basic>(*cit).tinfo()==this->tinfo()) {
			++nexpairseqs;
			noperands += ex_to<expairseq>(*cit).seq.size();
		}
		if (is_a<mul>(*this) && (!do_idx_rename) &&
				cit->info(info_flags::has_indices))
			do_idx_rename = true;
		++cit;
	}
	
	// reserve seq and coeffseq which will hold all operands
	seq.reserve(v.size()+noperands-nexpairseqs);
	
	// copy elements and split off numerical part
	make_flat_inserter mf(v, do_idx_rename);
	cit = v.begin();
	while (cit!=v.end()) {
		if (ex_to<basic>(*cit).tinfo()==this->tinfo()) {
			ex newfactor = mf.handle_factor(*cit, _ex1);
			const expairseq &subseqref = ex_to<expairseq>(newfactor);
			combine_overall_coeff(subseqref.overall_coeff);
			epvector::const_iterator cit_s = subseqref.seq.begin();
			while (cit_s!=subseqref.seq.end()) {
				seq.push_back(*cit_s);
				++cit_s;
			}
		} else {
			if (is_exactly_a<numeric>(*cit))
				combine_overall_coeff(*cit);
			else {
				ex newfactor = mf.handle_factor(*cit, _ex1);
				seq.push_back(split_ex_to_pair(newfactor));
			}
		}
		++cit;
	}
}

/** Combine this expairseq with argument epvector.
 *  It cares for associativity as well as for special handling of numerics. */
void expairseq::make_flat(const epvector &v, bool do_index_renaming)
{
	epvector::const_iterator cit;
	
	// count number of operands which are of same expairseq derived type
	// and their cumulative number of operands
	int nexpairseqs = 0;
	int noperands = 0;
	bool really_need_rename_inds = false;
	
	cit = v.begin();
	while (cit!=v.end()) {
		if (ex_to<basic>(cit->rest).tinfo()==this->tinfo()) {
			++nexpairseqs;
			noperands += ex_to<expairseq>(cit->rest).seq.size();
		}
		if ((!really_need_rename_inds) && is_a<mul>(*this) &&
				cit->rest.info(info_flags::has_indices))
			really_need_rename_inds = true;
		++cit;
	}
	do_index_renaming = do_index_renaming && really_need_rename_inds;
	
	// reserve seq and coeffseq which will hold all operands
	seq.reserve(v.size()+noperands-nexpairseqs);
	make_flat_inserter mf(v, do_index_renaming);
	
	// copy elements and split off numerical part
	cit = v.begin();
	while (cit!=v.end()) {
		if (ex_to<basic>(cit->rest).tinfo()==this->tinfo() &&
		    this->can_make_flat(*cit)) {
			ex newrest = mf.handle_factor(cit->rest, cit->coeff);
			const expairseq &subseqref = ex_to<expairseq>(newrest);
			combine_overall_coeff(ex_to<numeric>(subseqref.overall_coeff),
			                                    ex_to<numeric>(cit->coeff));
			epvector::const_iterator cit_s = subseqref.seq.begin();
			while (cit_s!=subseqref.seq.end()) {
				seq.push_back(expair(cit_s->rest,
				                     ex_to<numeric>(cit_s->coeff).mul_dyn(ex_to<numeric>(cit->coeff))));
				//seq.push_back(combine_pair_with_coeff_to_pair(*cit_s,
				//                                              (*cit).coeff));
				++cit_s;
			}
		} else {
			if (cit->is_canonical_numeric())
				combine_overall_coeff(mf.handle_factor(cit->rest, _ex1));
			else {
				ex rest = cit->rest;
				ex newrest = mf.handle_factor(rest, cit->coeff);
				if (are_ex_trivially_equal(newrest, rest))
					seq.push_back(*cit);
				else
					seq.push_back(expair(newrest, cit->coeff));
			}
		}
		++cit;
	}
}

/** Brings this expairseq into a sorted (canonical) form. */
void expairseq::canonicalize()
{
	std::sort(seq.begin(), seq.end(), expair_rest_is_less());
}


/** Compact a presorted expairseq by combining all matching expairs to one
 *  each.  On an add object, this is responsible for 2*x+3*x+y -> 5*x+y, for
 *  instance. */
void expairseq::combine_same_terms_sorted_seq()
{
	if (seq.size()<2)
		return;

	bool needs_further_processing = false;

	epvector::iterator itin1 = seq.begin();
	epvector::iterator itin2 = itin1+1;
	epvector::iterator itout = itin1;
	epvector::iterator last = seq.end();
	// must_copy will be set to true the first time some combination is 
	// possible from then on the sequence has changed and must be compacted
	bool must_copy = false;
	while (itin2!=last) {
		if (itin1->rest.compare(itin2->rest)==0) {
			itin1->coeff = ex_to<numeric>(itin1->coeff).
			               add_dyn(ex_to<numeric>(itin2->coeff));
			if (expair_needs_further_processing(itin1))
				needs_further_processing = true;
			must_copy = true;
		} else {
			if (!ex_to<numeric>(itin1->coeff).is_zero()) {
				if (must_copy)
					*itout = *itin1;
				++itout;
			}
			itin1 = itin2;
		}
		++itin2;
	}
	if (!ex_to<numeric>(itin1->coeff).is_zero()) {
		if (must_copy)
			*itout = *itin1;
		++itout;
	}
	if (itout!=last)
		seq.erase(itout,last);

	if (needs_further_processing) {
		epvector v = seq;
		seq.clear();
		construct_from_epvector(v);
	}
}

#if EXPAIRSEQ_USE_HASHTAB

unsigned expairseq::calc_hashtabsize(unsigned sz) const
{
	unsigned size;
	unsigned nearest_power_of_2 = 1 << log2(sz);
	// if (nearest_power_of_2 < maxhashtabsize/hashtabfactor) {
	//  size = nearest_power_of_2*hashtabfactor;
	size = nearest_power_of_2/hashtabfactor;
	if (size<minhashtabsize)
		return 0;

	// hashtabsize must be a power of 2
	GINAC_ASSERT((1U << log2(size))==size);
	return size;
}

unsigned expairseq::calc_hashindex(const ex &e) const
{
	// calculate hashindex
	unsigned hashindex;
	if (is_a<numeric>(e)) {
		hashindex = hashmask;
	} else {
		hashindex = e.gethash() & hashmask;
		// last hashtab entry is reserved for numerics
		if (hashindex==hashmask) hashindex = 0;
	}
	GINAC_ASSERT((hashindex<hashtabsize)||(hashtabsize==0));
	return hashindex;
}

void expairseq::shrink_hashtab()
{
	unsigned new_hashtabsize;
	while (hashtabsize!=(new_hashtabsize=calc_hashtabsize(seq.size()))) {
		GINAC_ASSERT(new_hashtabsize<hashtabsize);
		if (new_hashtabsize==0) {
			hashtab.clear();
			hashtabsize = 0;
			canonicalize();
			return;
		}
		
		// shrink by a factor of 2
		unsigned half_hashtabsize = hashtabsize/2;
		for (unsigned i=0; i<half_hashtabsize-1; ++i)
			hashtab[i].merge(hashtab[i+half_hashtabsize],epp_is_less());
		// special treatment for numeric hashes
		hashtab[0].merge(hashtab[half_hashtabsize-1],epp_is_less());
		hashtab[half_hashtabsize-1] = hashtab[hashtabsize-1];
		hashtab.resize(half_hashtabsize);
		hashtabsize = half_hashtabsize;
		hashmask = hashtabsize-1;
	}
}

void expairseq::remove_hashtab_entry(epvector::const_iterator element)
{
	if (hashtabsize==0)
		return; // nothing to do
	
	// calculate hashindex of element to be deleted
	unsigned hashindex = calc_hashindex((*element).rest);

	// find it in hashtab and remove it
	epplist &eppl = hashtab[hashindex];
	epplist::iterator epplit = eppl.begin();
	bool erased = false;
	while (epplit!=eppl.end()) {
		if (*epplit == element) {
			eppl.erase(epplit);
			erased = true;
			break;
		}
		++epplit;
	}
	if (!erased) {
		std::cout << "tried to erase " << element-seq.begin() << std::endl;
		std::cout << "size " << seq.end()-seq.begin() << std::endl;

		unsigned hashindex = calc_hashindex(element->rest);
		epplist &eppl = hashtab[hashindex];
		epplist::iterator epplit = eppl.begin();
		bool erased = false;
		while (epplit!=eppl.end()) {
			if (*epplit == element) {
				eppl.erase(epplit);
				erased = true;
				break;
			}
			++epplit;
		}
		GINAC_ASSERT(erased);
	}
	GINAC_ASSERT(erased);
}

void expairseq::move_hashtab_entry(epvector::const_iterator oldpos,
                                   epvector::iterator newpos)
{
	GINAC_ASSERT(hashtabsize!=0);
	
	// calculate hashindex of element which was moved
	unsigned hashindex=calc_hashindex((*newpos).rest);

	// find it in hashtab and modify it
	epplist &eppl = hashtab[hashindex];
	epplist::iterator epplit = eppl.begin();
	while (epplit!=eppl.end()) {
		if (*epplit == oldpos) {
			*epplit = newpos;
			break;
		}
		++epplit;
	}
	GINAC_ASSERT(epplit!=eppl.end());
}

void expairseq::sorted_insert(epplist &eppl, epvector::const_iterator elem)
{
	epplist::const_iterator current = eppl.begin();
	while ((current!=eppl.end()) && ((*current)->is_less(*elem))) {
		++current;
	}
	eppl.insert(current,elem);
}    

void expairseq::build_hashtab_and_combine(epvector::iterator &first_numeric,
                                          epvector::iterator &last_non_zero,
                                          std::vector<bool> &touched,
                                          unsigned &number_of_zeroes)
{
	epp current = seq.begin();

	while (current!=first_numeric) {
		if (is_exactly_a<numeric>(current->rest)) {
			--first_numeric;
			iter_swap(current,first_numeric);
		} else {
			// calculate hashindex
			unsigned currenthashindex = calc_hashindex(current->rest);

			// test if there is already a matching expair in the hashtab-list
			epplist &eppl=hashtab[currenthashindex];
			epplist::iterator epplit = eppl.begin();
			while (epplit!=eppl.end()) {
				if (current->rest.is_equal((*epplit)->rest))
					break;
				++epplit;
			}
			if (epplit==eppl.end()) {
				// no matching expair found, append this to end of list
				sorted_insert(eppl,current);
				++current;
			} else {
				// epplit points to a matching expair, combine it with current
				(*epplit)->coeff = ex_to<numeric>((*epplit)->coeff).
				                   add_dyn(ex_to<numeric>(current->coeff));
				
				// move obsolete current expair to end by swapping with last_non_zero element
				// if this was a numeric, it is swapped with the expair before first_numeric 
				iter_swap(current,last_non_zero);
				--first_numeric;
				if (first_numeric!=last_non_zero) iter_swap(first_numeric,current);
				--last_non_zero;
				++number_of_zeroes;
				// test if combined term has coeff 0 and can be removed is done later
				touched[(*epplit)-seq.begin()] = true;
			}
		}
	}
}

void expairseq::drop_coeff_0_terms(epvector::iterator &first_numeric,
                                   epvector::iterator &last_non_zero,
                                   std::vector<bool> &touched,
                                   unsigned &number_of_zeroes)
{
	// move terms with coeff 0 to end and remove them from hashtab
	// check only those elements which have been touched
	epp current = seq.begin();
	size_t i = 0;
	while (current!=first_numeric) {
		if (!touched[i]) {
			++current;
			++i;
		} else if (!ex_to<numeric>((*current).coeff).is_zero()) {
			++current;
			++i;
		} else {
			remove_hashtab_entry(current);
			
			// move element to the end, unless it is already at the end
			if (current!=last_non_zero) {
				iter_swap(current,last_non_zero);
				--first_numeric;
				bool numeric_swapped = first_numeric!=last_non_zero;
				if (numeric_swapped)
					iter_swap(first_numeric,current);
				epvector::iterator changed_entry;

				if (numeric_swapped)
					changed_entry = first_numeric;
				else
					changed_entry = last_non_zero;
				
				--last_non_zero;
				++number_of_zeroes;

				if (first_numeric!=current) {

					// change entry in hashtab which referred to first_numeric or last_non_zero to current
					move_hashtab_entry(changed_entry,current);
					touched[current-seq.begin()] = touched[changed_entry-seq.begin()];
				}
			} else {
				--first_numeric;
				--last_non_zero;
				++number_of_zeroes;
			}
		}
	}
	GINAC_ASSERT(i==current-seq.begin());
}

/** True if one of the coeffs vanishes, otherwise false.
 *  This would be an invariant violation, so this should only be used for
 *  debugging purposes. */
bool expairseq::has_coeff_0() const
{
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		if (i->coeff.is_zero())
			return true;
		++i;
	}
	return false;
}

void expairseq::add_numerics_to_hashtab(epvector::iterator first_numeric,
										epvector::const_iterator last_non_zero)
{
	if (first_numeric == seq.end()) return; // no numerics
	
	epvector::const_iterator current = first_numeric, last = last_non_zero + 1;
	while (current != last) {
		sorted_insert(hashtab[hashmask], current);
		++current;
	}
}

void expairseq::combine_same_terms()
{
	// combine same terms, drop term with coeff 0, move numerics to end
	
	// calculate size of hashtab
	hashtabsize = calc_hashtabsize(seq.size());
	
	// hashtabsize is a power of 2
	hashmask = hashtabsize-1;
	
	// allocate hashtab
	hashtab.clear();
	hashtab.resize(hashtabsize);
	
	if (hashtabsize==0) {
		canonicalize();
		combine_same_terms_sorted_seq();
		GINAC_ASSERT(!has_coeff_0());
		return;
	}
	
	// iterate through seq, move numerics to end,
	// fill hashtab and combine same terms
	epvector::iterator first_numeric = seq.end();
	epvector::iterator last_non_zero = seq.end()-1;
	
	size_t num = seq.size();
	std::vector<bool> touched(num);
	
	unsigned number_of_zeroes = 0;
	
	GINAC_ASSERT(!has_coeff_0());
	build_hashtab_and_combine(first_numeric,last_non_zero,touched,number_of_zeroes);
	
	// there should not be any terms with coeff 0 from the beginning,
	// so it should be safe to skip this step
	if (number_of_zeroes!=0) {
		drop_coeff_0_terms(first_numeric,last_non_zero,touched,number_of_zeroes);
	}
	
	add_numerics_to_hashtab(first_numeric,last_non_zero);
	
	// pop zero elements
	for (unsigned i=0; i<number_of_zeroes; ++i) {
		seq.pop_back();
	}
	
	// shrink hashtabsize to calculated value
	GINAC_ASSERT(!has_coeff_0());
	
	shrink_hashtab();
	
	GINAC_ASSERT(!has_coeff_0());
}

#endif // EXPAIRSEQ_USE_HASHTAB

/** Check if this expairseq is in sorted (canonical) form.  Useful mainly for
 *  debugging or in assertions since being sorted is an invariance. */
bool expairseq::is_canonical() const
{
	if (seq.size() <= 1)
		return 1;
	
#if EXPAIRSEQ_USE_HASHTAB
	if (hashtabsize > 0) return 1; // not canoncalized
#endif // EXPAIRSEQ_USE_HASHTAB
	
	epvector::const_iterator it = seq.begin(), itend = seq.end();
	epvector::const_iterator it_last = it;
	for (++it; it!=itend; it_last=it, ++it) {
		if (!(it_last->is_less(*it) || it_last->is_equal(*it))) {
			if (!is_exactly_a<numeric>(it_last->rest) ||
				!is_exactly_a<numeric>(it->rest)) {
				// double test makes it easier to set a breakpoint...
				if (!is_exactly_a<numeric>(it_last->rest) ||
					!is_exactly_a<numeric>(it->rest)) {
					printpair(std::clog, *it_last, 0);
					std::clog << ">";
					printpair(std::clog, *it, 0);
					std::clog << "\n";
					std::clog << "pair1:" << std::endl;
					it_last->rest.print(print_tree(std::clog));
					it_last->coeff.print(print_tree(std::clog));
					std::clog << "pair2:" << std::endl;
					it->rest.print(print_tree(std::clog));
					it->coeff.print(print_tree(std::clog));
					return 0;
				}
			}
		}
	}
	return 1;
}


/** Member-wise expand the expairs in this sequence.
 *
 *  @see expairseq::expand()
 *  @return pointer to epvector containing expanded pairs or zero pointer,
 *  if no members were changed. */
std::auto_ptr<epvector> expairseq::expandchildren(unsigned options) const
{
	const epvector::const_iterator last = seq.end();
	epvector::const_iterator cit = seq.begin();
	while (cit!=last) {
		const ex &expanded_ex = cit->rest.expand(options);
		if (!are_ex_trivially_equal(cit->rest,expanded_ex)) {
			
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
			s->push_back(combine_ex_with_coeff_to_pair(expanded_ex,
			                                           cit2->coeff));
			++cit2;

			// copy rest
			while (cit2!=last) {
				s->push_back(combine_ex_with_coeff_to_pair(cit2->rest.expand(options),
				                                           cit2->coeff));
				++cit2;
			}
			return s;
		}
		++cit;
	}
	
	return std::auto_ptr<epvector>(0); // signalling nothing has changed
}


/** Member-wise evaluate the expairs in this sequence.
 *
 *  @see expairseq::eval()
 *  @return pointer to epvector containing evaluated pairs or zero pointer,
 *  if no members were changed. */
std::auto_ptr<epvector> expairseq::evalchildren(int level) const
{
	// returns a NULL pointer if nothing had to be evaluated
	// returns a pointer to a newly created epvector otherwise
	// (which has to be deleted somewhere else)

	if (level==1)
		return std::auto_ptr<epvector>(0);
	
	if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	--level;
	epvector::const_iterator last = seq.end();
	epvector::const_iterator cit = seq.begin();
	while (cit!=last) {
		const ex &evaled_ex = cit->rest.eval(level);
		if (!are_ex_trivially_equal(cit->rest,evaled_ex)) {
			
			// something changed, copy seq, eval and return it
			std::auto_ptr<epvector> s(new epvector);
			s->reserve(seq.size());
			
			// copy parts of seq which are known not to have changed
			epvector::const_iterator cit2=seq.begin();
			while (cit2!=cit) {
				s->push_back(*cit2);
				++cit2;
			}

			// copy first changed element
			s->push_back(combine_ex_with_coeff_to_pair(evaled_ex,
			                                           cit2->coeff));
			++cit2;

			// copy rest
			while (cit2!=last) {
				s->push_back(combine_ex_with_coeff_to_pair(cit2->rest.eval(level),
				                                           cit2->coeff));
				++cit2;
			}
			return s;
		}
		++cit;
	}
	
	return std::auto_ptr<epvector>(0); // signalling nothing has changed
}

/** Member-wise substitute in this sequence.
 *
 *  @see expairseq::subs()
 *  @return pointer to epvector containing pairs after application of subs,
 *    or NULL pointer if no members were changed. */
std::auto_ptr<epvector> expairseq::subschildren(const exmap & m, unsigned options) const
{
	// When any of the objects to be substituted is a product or power
	// we have to recombine the pairs because the numeric coefficients may
	// be part of the search pattern.
	if (!(options & (subs_options::pattern_is_product | subs_options::pattern_is_not_product))) {

		// Search the list of substitutions and cache our findings
		for (exmap::const_iterator it = m.begin(); it != m.end(); ++it) {
			if (is_exactly_a<mul>(it->first) || is_exactly_a<power>(it->first)) {
				options |= subs_options::pattern_is_product;
				break;
			}
		}
		if (!(options & subs_options::pattern_is_product))
			options |= subs_options::pattern_is_not_product;
	}

	if (options & subs_options::pattern_is_product) {

		// Substitute in the recombined pairs
		epvector::const_iterator cit = seq.begin(), last = seq.end();
		while (cit != last) {

			const ex &orig_ex = recombine_pair_to_ex(*cit);
			const ex &subsed_ex = orig_ex.subs(m, options);
			if (!are_ex_trivially_equal(orig_ex, subsed_ex)) {

				// Something changed, copy seq, subs and return it
				std::auto_ptr<epvector> s(new epvector);
				s->reserve(seq.size());

				// Copy parts of seq which are known not to have changed
				s->insert(s->begin(), seq.begin(), cit);

				// Copy first changed element
				s->push_back(split_ex_to_pair(subsed_ex));
				++cit;

				// Copy rest
				while (cit != last) {
					s->push_back(split_ex_to_pair(recombine_pair_to_ex(*cit).subs(m, options)));
					++cit;
				}
				return s;
			}

			++cit;
		}

	} else {

		// Substitute only in the "rest" part of the pairs
		epvector::const_iterator cit = seq.begin(), last = seq.end();
		while (cit != last) {

			const ex &subsed_ex = cit->rest.subs(m, options);
			if (!are_ex_trivially_equal(cit->rest, subsed_ex)) {
			
				// Something changed, copy seq, subs and return it
				std::auto_ptr<epvector> s(new epvector);
				s->reserve(seq.size());

				// Copy parts of seq which are known not to have changed
				s->insert(s->begin(), seq.begin(), cit);
			
				// Copy first changed element
				s->push_back(combine_ex_with_coeff_to_pair(subsed_ex, cit->coeff));
				++cit;

				// Copy rest
				while (cit != last) {
					s->push_back(combine_ex_with_coeff_to_pair(cit->rest.subs(m, options), cit->coeff));
					++cit;
				}
				return s;
			}

			++cit;
		}
	}
	
	// Nothing has changed
	return std::auto_ptr<epvector>(0);
}

//////////
// static member variables
//////////

#if EXPAIRSEQ_USE_HASHTAB
unsigned expairseq::maxhashtabsize = 0x4000000U;
unsigned expairseq::minhashtabsize = 0x1000U;
unsigned expairseq::hashtabfactor = 1;
#endif // EXPAIRSEQ_USE_HASHTAB

} // namespace GiNaC
