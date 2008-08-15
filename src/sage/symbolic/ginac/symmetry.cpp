/** @file symmetry.cpp
 *
 *  Implementation of GiNaC's symmetry definitions. */

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
#include <stdexcept>
#include <functional>
#include <limits>

#include "symmetry.h"
#include "lst.h"
#include "numeric.h" // for factorial()
#include "operators.h"
#include "archive.h"
#include "utils.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(symmetry, basic,
  print_func<print_context>(&symmetry::do_print).
  print_func<print_tree>(&symmetry::do_print_tree))

/*
   Some notes about the structure of a symmetry tree:
    - The leaf nodes of the tree are of type "none", have one index, and no
      children (of course). They are constructed by the symmetry(unsigned)
      constructor.
    - Leaf nodes are the only nodes that only have one index.
    - Container nodes contain two or more children. The "indices" set member
      is the set union of the index sets of all children, and the "children"
      vector stores the children themselves.
    - The index set of each child of a "symm", "anti" or "cycl" node must
      have the same size. It follows that the children of such a node are
      either all leaf nodes, or all container nodes with two or more indices.
*/

//////////
// default constructor
//////////

symmetry::symmetry() : inherited(&symmetry::tinfo_static), type(none)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

symmetry::symmetry(unsigned i) : inherited(&symmetry::tinfo_static), type(none)
{
	indices.insert(i);
	setflag(status_flags::evaluated | status_flags::expanded);
}

symmetry::symmetry(symmetry_type t, const symmetry &c1, const symmetry &c2) : inherited(&symmetry::tinfo_static), type(t)
{
	add(c1); add(c2);
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// archiving
//////////

/** Construct object from archive_node. */
symmetry::symmetry(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	unsigned t;
	if (!(n.find_unsigned("type", t)))
		throw (std::runtime_error("unknown symmetry type in archive"));
	type = (symmetry_type)t;

	unsigned i = 0;
	while (true) {
		ex e;
		if (n.find_ex("child", e, sym_lst, i))
			add(ex_to<symmetry>(e));
		else
			break;
		i++;
	}

	if (i == 0) {
		while (true) {
			unsigned u;
			if (n.find_unsigned("index", u, i))
				indices.insert(u);
			else
				break;
			i++;
		}
	}
}

/** Archive the object. */
void symmetry::archive(archive_node &n) const
{
	inherited::archive(n);

	n.add_unsigned("type", type);

	if (children.empty()) {
		std::set<unsigned>::const_iterator i = indices.begin(), iend = indices.end();
		while (i != iend) {
			n.add_unsigned("index", *i);
			i++;
		}
	} else {
		exvector::const_iterator i = children.begin(), iend = children.end();
		while (i != iend) {
			n.add_ex("child", *i);
			i++;
		}
	}
}

DEFAULT_UNARCHIVE(symmetry)

//////////
// functions overriding virtual functions from base classes
//////////

int symmetry::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<symmetry>(other));

	// For archiving purposes we need to have an ordering of symmetries.
	const symmetry &othersymm = ex_to<symmetry>(other);

	// Compare type.
	if (type > othersymm.type)
		return 1;
	if (type < othersymm.type)
		return -1;

	// Compare the index set.
	size_t this_size = indices.size();
	size_t that_size = othersymm.indices.size();
	if (this_size > that_size)
		return 1;
	if (this_size < that_size)
		return -1;
	typedef std::set<unsigned>::iterator set_it;
	set_it end = indices.end();
	for (set_it i=indices.begin(),j=othersymm.indices.begin(); i!=end; ++i,++j) {
		if(*i < *j)
			return 1;
		if(*i > *j)
			return -1;
	}

	// Compare the children.
	if (children.size() > othersymm.children.size())
		return 1;
	if (children.size() < othersymm.children.size())
		return -1;
	for (size_t i=0; i<children.size(); ++i) {
		int cmpval = ex_to<symmetry>(children[i])
			.compare_same_type(ex_to<symmetry>(othersymm.children[i]));
		if (cmpval)
			return cmpval;
	}

	return 0;
}

unsigned symmetry::calchash() const
{
	unsigned v = golden_ratio_hash((p_int)tinfo());

	if (type == none) {
		v = rotate_left(v);
		v ^= *(indices.begin());
	} else {
		for (exvector::const_iterator i=children.begin(); i!=children.end(); ++i)
		{
			v = rotate_left(v);
			v ^= i->gethash();
		}
	}

	if (flags & status_flags::evaluated) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}

	return v;
}

void symmetry::do_print(const print_context & c, unsigned level) const
{
	if (children.empty()) {
		if (indices.size() > 0)
			c.s << *(indices.begin());
		else
			c.s << "none";
	} else {
		switch (type) {
			case none: c.s << '!'; break;
			case symmetric: c.s << '+'; break;
			case antisymmetric: c.s << '-'; break;
			case cyclic: c.s << '@'; break;
			default: c.s << '?'; break;
		}
		c.s << '(';
		size_t num = children.size();
		for (size_t i=0; i<num; i++) {
			children[i].print(c);
			if (i != num - 1)
				c.s << ",";
		}
		c.s << ')';
	}
}

void symmetry::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", type=";

	switch (type) {
		case none: c.s << "none"; break;
		case symmetric: c.s << "symm"; break;
		case antisymmetric: c.s << "anti"; break;
		case cyclic: c.s << "cycl"; break;
		default: c.s << "<unknown>"; break;
	}

	c.s << ", indices=(";
	if (!indices.empty()) {
		std::set<unsigned>::const_iterator i = indices.begin(), end = indices.end();
		--end;
		while (i != end)
			c.s << *i++ << ",";
		c.s << *i;
	}
	c.s << ")\n";

	exvector::const_iterator i = children.begin(), end = children.end();
	while (i != end) {
		i->print(c, level + c.delta_indent);
		++i;
	}
}

//////////
// non-virtual functions in this class
//////////

bool symmetry::has_cyclic() const
{
	if (type == cyclic)
		return true;

	for (exvector::const_iterator i=children.begin(); i!=children.end(); ++i)
		if (ex_to<symmetry>(*i).has_cyclic())
			return true;

	return false;
}

symmetry &symmetry::add(const symmetry &c)
{
	// All children must have the same number of indices
	if (type != none && !children.empty()) {
		GINAC_ASSERT(is_exactly_a<symmetry>(children[0]));
		if (ex_to<symmetry>(children[0]).indices.size() != c.indices.size())
			throw (std::logic_error("symmetry:add(): children must have same number of indices"));
	}

	// Compute union of indices and check whether the two sets are disjoint
	std::set<unsigned> un;
	set_union(indices.begin(), indices.end(), c.indices.begin(), c.indices.end(), inserter(un, un.begin()));
	if (un.size() != indices.size() + c.indices.size())
		throw (std::logic_error("symmetry::add(): the same index appears in more than one child"));

	// Set new index set
	indices.swap(un);

	// Add child node
	children.push_back(c);
	return *this;
}

void symmetry::validate(unsigned n)
{
	if (indices.upper_bound(n - 1) != indices.end())
		throw (std::range_error("symmetry::verify(): index values are out of range"));
	if (type != none && indices.empty()) {
		for (unsigned i=0; i<n; i++)
			add(i);
	}
}

//////////
// global functions
//////////

static const symmetry & index0()
{
	static ex s = (new symmetry(0))->setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

static const symmetry & index1()
{
	static ex s = (new symmetry(1))->setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

static const symmetry & index2()
{
	static ex s = (new symmetry(2))->setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

static const symmetry & index3()
{
	static ex s = (new symmetry(3))->setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

const symmetry & not_symmetric()
{
	static ex s = (new symmetry)->setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

const symmetry & symmetric2()
{
	static ex s = (new symmetry(symmetry::symmetric, index0(), index1()))->setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

const symmetry & symmetric3()
{
	static ex s = (new symmetry(symmetry::symmetric, index0(), index1()))->add(index2()).setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

const symmetry & symmetric4()
{
	static ex s = (new symmetry(symmetry::symmetric, index0(), index1()))->add(index2()).add(index3()).setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

const symmetry & antisymmetric2()
{
	static ex s = (new symmetry(symmetry::antisymmetric, index0(), index1()))->setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

const symmetry & antisymmetric3()
{
	static ex s = (new symmetry(symmetry::antisymmetric, index0(), index1()))->add(index2()).setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

const symmetry & antisymmetric4()
{
	static ex s = (new symmetry(symmetry::antisymmetric, index0(), index1()))->add(index2()).add(index3()).setflag(status_flags::dynallocated);
	return ex_to<symmetry>(s);
}

class sy_is_less : public std::binary_function<ex, ex, bool> {
	exvector::iterator v;

public:
	sy_is_less(exvector::iterator v_) : v(v_) {}

	bool operator() (const ex &lh, const ex &rh) const
	{
		GINAC_ASSERT(is_exactly_a<symmetry>(lh));
		GINAC_ASSERT(is_exactly_a<symmetry>(rh));
		GINAC_ASSERT(ex_to<symmetry>(lh).indices.size() == ex_to<symmetry>(rh).indices.size());
		std::set<unsigned>::const_iterator ait = ex_to<symmetry>(lh).indices.begin(), aitend = ex_to<symmetry>(lh).indices.end(), bit = ex_to<symmetry>(rh).indices.begin();
		while (ait != aitend) {
			int cmpval = v[*ait].compare(v[*bit]);
			if (cmpval < 0)
				return true;
			else if (cmpval > 0)
				return false;
			++ait; ++bit;
		}
		return false;
	}
};

class sy_swap : public std::binary_function<ex, ex, void> {
	exvector::iterator v;

public:
	bool &swapped;

	sy_swap(exvector::iterator v_, bool &s) : v(v_), swapped(s) {}

	void operator() (const ex &lh, const ex &rh)
	{
		GINAC_ASSERT(is_exactly_a<symmetry>(lh));
		GINAC_ASSERT(is_exactly_a<symmetry>(rh));
		GINAC_ASSERT(ex_to<symmetry>(lh).indices.size() == ex_to<symmetry>(rh).indices.size());
		std::set<unsigned>::const_iterator ait = ex_to<symmetry>(lh).indices.begin(), aitend = ex_to<symmetry>(lh).indices.end(), bit = ex_to<symmetry>(rh).indices.begin();
		while (ait != aitend) {
			v[*ait].swap(v[*bit]);
			++ait; ++bit;
		}
		swapped = true;
	}
};

int canonicalize(exvector::iterator v, const symmetry &symm)
{
	// Less than two elements? Then do nothing
	if (symm.indices.size() < 2)
		return std::numeric_limits<int>::max();

	// Canonicalize children first
	bool something_changed = false;
	int sign = 1;
	exvector::const_iterator first = symm.children.begin(), last = symm.children.end();
	while (first != last) {
		GINAC_ASSERT(is_exactly_a<symmetry>(*first));
		int child_sign = canonicalize(v, ex_to<symmetry>(*first));
		if (child_sign == 0)
			return 0;
		if (child_sign != std::numeric_limits<int>::max()) {
			something_changed = true;
			sign *= child_sign;
		}
		first++;
	}

	// Now reorder the children
	first = symm.children.begin();
	switch (symm.type) {
		case symmetry::symmetric:
			// Sort the children in ascending order
			shaker_sort(first, last, sy_is_less(v), sy_swap(v, something_changed));
			break;
		case symmetry::antisymmetric:
			// Sort the children in ascending order, keeping track of the signum
			sign *= permutation_sign(first, last, sy_is_less(v), sy_swap(v, something_changed));
			if (sign == 0)
				return 0;
			break;
		case symmetry::cyclic:
			// Permute the smallest child to the front
			cyclic_permutation(first, last, min_element(first, last, sy_is_less(v)), sy_swap(v, something_changed));
			break;
		default:
			break;
	}
	return something_changed ? sign : std::numeric_limits<int>::max();
}


// Symmetrize/antisymmetrize over a vector of objects
static ex symm(const ex & e, exvector::const_iterator first, exvector::const_iterator last, bool asymmetric)
{
	// Need at least 2 objects for this operation
	unsigned num = last - first;
	if (num < 2)
		return e;

	// Transform object vector to a lst (for subs())
	lst orig_lst(first, last);

	// Create index vectors for permutation
	unsigned *iv = new unsigned[num], *iv2;
	for (unsigned i=0; i<num; i++)
		iv[i] = i;
	iv2 = (asymmetric ? new unsigned[num] : NULL);

	// Loop over all permutations (the first permutation, which is the
	// identity, is unrolled)
	ex sum = e;
	while (std::next_permutation(iv, iv + num)) {
		lst new_lst;
		for (unsigned i=0; i<num; i++)
			new_lst.append(orig_lst.op(iv[i]));
		ex term = e.subs(orig_lst, new_lst, subs_options::no_pattern|subs_options::no_index_renaming);
		if (asymmetric) {
			memcpy(iv2, iv, num * sizeof(unsigned));
			term *= permutation_sign(iv2, iv2 + num);
		}
		sum += term;
	}

	delete[] iv;
	delete[] iv2;

	return sum / factorial(numeric(num));
}

ex symmetrize(const ex & e, exvector::const_iterator first, exvector::const_iterator last)
{
	return symm(e, first, last, false);
}

ex antisymmetrize(const ex & e, exvector::const_iterator first, exvector::const_iterator last)
{
	return symm(e, first, last, true);
}

ex symmetrize_cyclic(const ex & e, exvector::const_iterator first, exvector::const_iterator last)
{
	// Need at least 2 objects for this operation
	unsigned num = last - first;
	if (num < 2)
		return e;

	// Transform object vector to a lst (for subs())
	lst orig_lst(first, last);
	lst new_lst = orig_lst;

	// Loop over all cyclic permutations (the first permutation, which is
	// the identity, is unrolled)
	ex sum = e;
	for (unsigned i=0; i<num-1; i++) {
		ex perm = new_lst.op(0);
		new_lst.remove_first().append(perm);
		sum += e.subs(orig_lst, new_lst, subs_options::no_pattern|subs_options::no_index_renaming);
	}
	return sum / num;
}

/** Symmetrize expression over a list of objects (symbols, indices). */
ex ex::symmetrize(const lst & l) const
{
	exvector v(l.begin(), l.end());
	return symm(*this, v.begin(), v.end(), false);
}

/** Antisymmetrize expression over a list of objects (symbols, indices). */
ex ex::antisymmetrize(const lst & l) const
{
	exvector v(l.begin(), l.end());
	return symm(*this, v.begin(), v.end(), true);
}

/** Symmetrize expression by cyclic permutation over a list of objects
 *  (symbols, indices). */
ex ex::symmetrize_cyclic(const lst & l) const
{
	exvector v(l.begin(), l.end());
	return GiNaC::symmetrize_cyclic(*this, v.begin(), v.end());
}

} // namespace GiNaC
