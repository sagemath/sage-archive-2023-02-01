/** @file idx.cpp
 *
 *  Implementation of GiNaC's indices. */

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
#include <sstream>
#include <stdexcept>

#include "idx.h"
#include "symbol.h"
#include "lst.h"
#include "relational.h"
#include "operators.h"
#include "archive.h"
#include "utils.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(idx, basic,
  print_func<print_context>(&idx::do_print).
  print_func<print_latex>(&idx::do_print_latex).
  print_func<print_csrc>(&idx::do_print_csrc).
  print_func<print_tree>(&idx::do_print_tree))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(varidx, idx,
  print_func<print_context>(&varidx::do_print).
  print_func<print_latex>(&varidx::do_print_latex).
  print_func<print_tree>(&varidx::do_print_tree))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(spinidx, varidx,
  print_func<print_context>(&spinidx::do_print).
  print_func<print_latex>(&spinidx::do_print_latex).
  print_func<print_tree>(&spinidx::do_print_tree))

//////////
// default constructor
//////////

idx::idx() : inherited(&idx::tinfo_static) {}

varidx::varidx() : covariant(false)
{
	tinfo_key = &varidx::tinfo_static;
}

spinidx::spinidx() : dotted(false)
{
	tinfo_key = &spinidx::tinfo_static;
}

//////////
// other constructors
//////////

idx::idx(const ex & v, const ex & d) : inherited(&idx::tinfo_static), value(v), dim(d)
{
	if (is_dim_numeric())
		if (!dim.info(info_flags::posint))
			throw(std::invalid_argument("dimension of space must be a positive integer"));
}

varidx::varidx(const ex & v, const ex & d, bool cov) : inherited(v, d), covariant(cov)
{
	tinfo_key = &varidx::tinfo_static;
}

spinidx::spinidx(const ex & v, const ex & d, bool cov, bool dot) : inherited(v, d, cov), dotted(dot)
{
	tinfo_key = &spinidx::tinfo_static;
}

//////////
// archiving
//////////

idx::idx(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	n.find_ex("value", value, sym_lst);
	n.find_ex("dim", dim, sym_lst);
}

varidx::varidx(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	n.find_bool("covariant", covariant);
}

spinidx::spinidx(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	n.find_bool("dotted", dotted);
}

void idx::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_ex("value", value);
	n.add_ex("dim", dim);
}

void varidx::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_bool("covariant", covariant);
}

void spinidx::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_bool("dotted", dotted);
}

DEFAULT_UNARCHIVE(idx)
DEFAULT_UNARCHIVE(varidx)
DEFAULT_UNARCHIVE(spinidx)

//////////
// functions overriding virtual functions from base classes
//////////

void idx::print_index(const print_context & c, unsigned level) const
{
	bool need_parens = !(is_exactly_a<numeric>(value) || is_a<symbol>(value));
	if (need_parens)
		c.s << "(";
	value.print(c);
	if (need_parens)
		c.s << ")";
	if (c.options & print_options::print_index_dimensions) {
		c.s << "[";
		dim.print(c);
		c.s << "]";
	}
}

void idx::do_print(const print_context & c, unsigned level) const
{
	c.s << ".";
	print_index(c, level);
}

void idx::do_print_latex(const print_latex & c, unsigned level) const
{
	c.s << "{";
	print_index(c, level);
	c.s << "}";
}

void idx::do_print_csrc(const print_csrc & c, unsigned level) const
{
	c.s << "[";
	if (value.info(info_flags::integer))
		c.s << ex_to<numeric>(value).to_int();
	else
		value.print(c);
	c.s << "]";
}

void idx::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << std::endl;
	value.print(c, level +  c.delta_indent);
	dim.print(c, level + c.delta_indent);
}

void varidx::do_print(const print_context & c, unsigned level) const
{
	if (covariant)
		c.s << ".";
	else
		c.s << "~";
	print_index(c, level);
}

void varidx::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << (covariant ? ", covariant" : ", contravariant")
	    << std::endl;
	value.print(c, level + c.delta_indent);
	dim.print(c, level + c.delta_indent);
}

void spinidx::do_print(const print_context & c, unsigned level) const
{
	if (covariant)
		c.s << ".";
	else
		c.s << "~";
	if (dotted)
		c.s << "*";
	print_index(c, level);
}

void spinidx::do_print_latex(const print_latex & c, unsigned level) const
{
	if (dotted)
		c.s << "\\dot{";
	else
		c.s << "{";
	print_index(c, level);
	c.s << "}";
}

void spinidx::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << (covariant ? ", covariant" : ", contravariant")
	    << (dotted ? ", dotted" : ", undotted")
	    << std::endl;
	value.print(c, level + c.delta_indent);
	dim.print(c, level + c.delta_indent);
}

bool idx::info(unsigned inf) const
{
	switch(inf) {
		case info_flags::idx:
		case info_flags::has_indices:
			return true;
	}
	return inherited::info(inf);
}

size_t idx::nops() const
{
	// don't count the dimension as that is not really a sub-expression
	return 1;
}

ex idx::op(size_t i) const
{
	GINAC_ASSERT(i == 0);
	return value;
}

ex idx::map(map_function & f) const
{
	const ex &mapped_value = f(value);
	if (are_ex_trivially_equal(value, mapped_value))
		return *this;
	else {
		idx *copy = duplicate();
		copy->setflag(status_flags::dynallocated);
		copy->clearflag(status_flags::hash_calculated);
		copy->value = mapped_value;
		return *copy;
	}
}

/** Returns order relation between two indices of the same type. The order
 *  must be such that dummy indices lie next to each other. */
int idx::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<idx>(other));
	const idx &o = static_cast<const idx &>(other);

	int cmpval = value.compare(o.value);
	if (cmpval)
		return cmpval;
	return dim.compare(o.dim);
}

bool idx::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<idx>(other));
	const idx &o = static_cast<const idx &>(other);

	return dim.is_equal(o.dim);
}

int varidx::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<varidx>(other));
	const varidx &o = static_cast<const varidx &>(other);

	int cmpval = inherited::compare_same_type(other);
	if (cmpval)
		return cmpval;

	// Check variance last so dummy indices will end up next to each other
	if (covariant != o.covariant)
		return covariant ? -1 : 1;

	return 0;
}

bool varidx::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<varidx>(other));
	const varidx &o = static_cast<const varidx &>(other);

	if (covariant != o.covariant)
		return false;

	return inherited::match_same_type(other);
}

int spinidx::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<spinidx>(other));
	const spinidx &o = static_cast<const spinidx &>(other);

	// Check dottedness first so dummy indices will end up next to each other
	if (dotted != o.dotted)
		return dotted ? -1 : 1;

	int cmpval = inherited::compare_same_type(other);
	if (cmpval)
		return cmpval;

	return 0;
}

bool spinidx::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<spinidx>(other));
	const spinidx &o = static_cast<const spinidx &>(other);

	if (dotted != o.dotted)
		return false;
	return inherited::match_same_type(other);
}

unsigned idx::calchash() const
{
	// NOTE: The code in simplify_indexed() assumes that canonically
	// ordered sequences of indices have the two members of dummy index
	// pairs lying next to each other. The hash values for indices must
	// be devised accordingly. The easiest (only?) way to guarantee the
	// desired ordering is to make indices with the same value have equal
	// hash keys. That is, the hash values must not depend on the index
	// dimensions or other attributes (variance etc.).
	// The compare_same_type() methods will take care of the rest.
	unsigned v = golden_ratio_hash((p_int)tinfo());
	v = rotate_left(v);
	v ^= value.gethash();

	// Store calculated hash value only if object is already evaluated
	if (flags & status_flags::evaluated) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}

	return v;
}

/** By default, basic::evalf would evaluate the index value but we don't want
 *  a.1 to become a.(1.0). */
ex idx::evalf(int level, int prec) const
{
	return *this;
}

ex idx::subs(const exmap & m, unsigned options) const
{
	// First look for index substitutions
	exmap::const_iterator it = m.find(*this);
	if (it != m.end()) {

		// Substitution index->index
		if (is_a<idx>(it->second) || (options & subs_options::really_subs_idx))
			return it->second;

		// Otherwise substitute value
		idx *i_copy = duplicate();
		i_copy->value = it->second;
		i_copy->clearflag(status_flags::hash_calculated);
		return i_copy->setflag(status_flags::dynallocated);
	}

	// None, substitute objects in value (not in dimension)
	const ex &subsed_value = value.subs(m, options);
	if (are_ex_trivially_equal(value, subsed_value))
		return *this;

	idx *i_copy = duplicate();
	i_copy->value = subsed_value;
	i_copy->clearflag(status_flags::hash_calculated);
	return i_copy->setflag(status_flags::dynallocated);
}

/** Implementation of ex::diff() for an index always returns 0.
 *
 *  @see ex::diff */
ex idx::derivative(const symbol & s) const
{
	return _ex0;
}

//////////
// new virtual functions
//////////

bool idx::is_dummy_pair_same_type(const basic & other) const
{
	const idx &o = static_cast<const idx &>(other);

	// Only pure symbols form dummy pairs, "2n+1" doesn't
	if (!is_a<symbol>(value))
		return false;

	// Value must be equal, of course
	if (!value.is_equal(o.value))
		return false;

	// Dimensions need not be equal but must be comparable (so we can
	// determine the minimum dimension of contractions)
	if (dim.is_equal(o.dim))
		return true;

	return is_exactly_a<numeric>(dim) || is_exactly_a<numeric>(o.dim);
}

bool varidx::is_dummy_pair_same_type(const basic & other) const
{
	const varidx &o = static_cast<const varidx &>(other);

	// Variance must be opposite
	if (covariant == o.covariant)
		return false;

	return inherited::is_dummy_pair_same_type(other);
}

bool spinidx::is_dummy_pair_same_type(const basic & other) const
{
	const spinidx &o = static_cast<const spinidx &>(other);

	// Dottedness must be the same
	if (dotted != o.dotted)
		return false;

	return inherited::is_dummy_pair_same_type(other);
}


//////////
// non-virtual functions
//////////

ex idx::replace_dim(const ex & new_dim) const
{
	idx *i_copy = duplicate();
	i_copy->dim = new_dim;
	i_copy->clearflag(status_flags::hash_calculated);
	return i_copy->setflag(status_flags::dynallocated);
}

ex idx::minimal_dim(const idx & other) const
{
	return GiNaC::minimal_dim(dim, other.dim);
}

ex varidx::toggle_variance() const
{
	varidx *i_copy = duplicate();
	i_copy->covariant = !i_copy->covariant;
	i_copy->clearflag(status_flags::hash_calculated);
	return i_copy->setflag(status_flags::dynallocated);
}

ex spinidx::toggle_dot() const
{
	spinidx *i_copy = duplicate();
	i_copy->dotted = !i_copy->dotted;
	i_copy->clearflag(status_flags::hash_calculated);
	return i_copy->setflag(status_flags::dynallocated);
}

ex spinidx::toggle_variance_dot() const
{
	spinidx *i_copy = duplicate();
	i_copy->covariant = !i_copy->covariant;
	i_copy->dotted = !i_copy->dotted;
	i_copy->clearflag(status_flags::hash_calculated);
	return i_copy->setflag(status_flags::dynallocated);
}

//////////
// global functions
//////////

bool is_dummy_pair(const idx & i1, const idx & i2)
{
	// The indices must be of exactly the same type
	if (i1.tinfo() != i2.tinfo())
		return false;

	// Same type, let the indices decide whether they are paired
	return i1.is_dummy_pair_same_type(i2);
}

bool is_dummy_pair(const ex & e1, const ex & e2)
{
	// The expressions must be indices
	if (!is_a<idx>(e1) || !is_a<idx>(e2))
		return false;

	return is_dummy_pair(ex_to<idx>(e1), ex_to<idx>(e2));
}

void find_free_and_dummy(exvector::const_iterator it, exvector::const_iterator itend, exvector & out_free, exvector & out_dummy)
{
	out_free.clear();
	out_dummy.clear();

	// No indices? Then do nothing
	if (it == itend)
		return;

	// Only one index? Then it is a free one if it's not numeric
	if (itend - it == 1) {
		if (ex_to<idx>(*it).is_symbolic())
			out_free.push_back(*it);
		return;
	}

	// Sort index vector. This will cause dummy indices come to lie next
	// to each other (because the sort order is defined to guarantee this).
	exvector v(it, itend);
	shaker_sort(v.begin(), v.end(), ex_is_less(), ex_swap());

	// Find dummy pairs and free indices
	it = v.begin(); itend = v.end();
	exvector::const_iterator last = it++;
	while (it != itend) {
		if (is_dummy_pair(*it, *last)) {
			out_dummy.push_back(*last);
			it++;
			if (it == itend)
				return;
		} else {
			if (!it->is_equal(*last) && ex_to<idx>(*last).is_symbolic())
				out_free.push_back(*last);
		}
		last = it++;
	}
	if (ex_to<idx>(*last).is_symbolic())
		out_free.push_back(*last);
}

ex minimal_dim(const ex & dim1, const ex & dim2)
{
	if (dim1.is_equal(dim2) || dim1 < dim2 || (is_exactly_a<numeric>(dim1) && !is_a<numeric>(dim2)))
		return dim1;
	else if (dim1 > dim2 || (!is_a<numeric>(dim1) && is_exactly_a<numeric>(dim2)))
		return dim2;
	else {
		std::ostringstream s;
		s << "minimal_dim(): index dimensions " << dim1 << " and " << dim2 << " cannot be ordered";
		throw (std::runtime_error(s.str()));
	}
}

} // namespace GiNaC
