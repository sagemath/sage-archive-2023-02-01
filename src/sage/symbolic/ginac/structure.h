/** @file structure.h
 *
 *  Wrapper template for making GiNaC classes out of C++ structures. */

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

#ifndef __GINAC_STRUCTURE_H__
#define __GINAC_STRUCTURE_H__

#include <functional>

#include "ex.h"
#include "ncmul.h"
#include "numeric.h"
#include "operators.h"
#include "print.h"

namespace GiNaC {


/** Comparison policy: all structures of one type are equal */
template <class T>
class compare_all_equal {
protected:
	static bool struct_is_equal(const T * t1, const T * t2) { return true; }
	static int struct_compare(const T * t1, const T * t2) { return 0; }

	// disallow destruction of structure through a compare_all_equal*
protected:
	~compare_all_equal() {}
};


/** Comparison policy: use std::equal_to/std::less (defaults to operators
 *  == and <) to compare structures. */
template <class T>
class compare_std_less {
protected:
	static bool struct_is_equal(const T * t1, const T * t2)
	{
		return std::equal_to<T>()(*t1, *t2);
	}

	static int struct_compare(const T * t1, const T * t2)
	{
		if (std::less<T>()(*t1, *t2))
			return -1;
		else if (std::less<T>()(*t2, *t1))
			return 1;
		else
			return 0;
	}

	// disallow destruction of structure through a compare_std_less*
protected:
	~compare_std_less() {}
};


/** Comparison policy: use bit-wise comparison to compare structures. */
template <class T>
class compare_bitwise {
protected:
	static bool struct_is_equal(const T * t1, const T * t2)
	{
		const char * cp1 = reinterpret_cast<const char *>(t1);
		const char * cp2 = reinterpret_cast<const char *>(t2);

		return std::equal(cp1, cp1 + sizeof(T), cp2);
	}

	static int struct_compare(const T * t1, const T * t2)
	{
		const unsigned char * cp1 = reinterpret_cast<const unsigned char *>(t1);
		const unsigned char * cp2 = reinterpret_cast<const unsigned char *>(t2);
		typedef std::pair<const unsigned char *, const unsigned char *> cppair;

		cppair res = std::mismatch(cp1, cp1 + sizeof(T), cp2);

		if (res.first == cp1 + sizeof(T))
			return 0;
		else if (*res.first < *res.second)
			return -1;
		else
			return 1;
	}

	// disallow destruction of structure through a compare_bitwise*
protected:
	~compare_bitwise() {}
};


// Select default comparison policy
template <class T, template <class> class ComparisonPolicy = compare_all_equal> class structure;


/** Wrapper template for making GiNaC classes out of C++ structures. */
template <class T, template <class> class ComparisonPolicy>
class structure : public basic, public ComparisonPolicy<T> {
	GINAC_DECLARE_REGISTERED_CLASS(structure, basic)

	// helpers
	static tinfo_t get_tinfo() { return reg_info.options.get_id(); }
	static const char *get_class_name() { return "structure"; }
	static tinfo_t next_structure_tinfo_key() { return new tinfo_static_t; }

	// constructors
public:
	/** Construct structure as a copy of a given C++ structure. */
	structure(const T & t) : inherited(get_tinfo()), obj(t) { }

	// functions overriding virtual functions from base classes
	// All these are just defaults that can be specialized by the user
public:
	// evaluation
	ex eval(int level = 0) const { return hold(); }
	ex evalf(int level = 0) const { return inherited::evalf(level); }
	ex evalm() const { return inherited::evalm(); }
protected:
	ex eval_ncmul(const exvector & v) const { return hold_ncmul(v); }
public:
	ex eval_indexed(const basic & i) const { return i.hold(); }

	// printing
	void print(const print_context & c, unsigned level = 0) const { inherited::print(c, level); }
	unsigned precedence() const { return 70; }

	// info
	bool info(unsigned inf) const { return false; }

	// operand access
	size_t nops() const { return 0; }
	ex op(size_t i) const { return inherited::op(i); }
	ex operator[](const ex & index) const { return inherited::operator[](index); }
	ex operator[](size_t i) const { return inherited::operator[](i); }
	ex & let_op(size_t i) { return inherited::let_op(i); }
	ex & operator[](const ex & index) { return inherited::operator[](index); }
	ex & operator[](size_t i) { return inherited::operator[](i); }

	// pattern matching
	bool has(const ex & other, unsigned options = 0) const { return inherited::has(other, options); }
	bool match(const ex & pattern, lst & repl_lst) const { return inherited::match(pattern, repl_lst); }
protected:
	bool match_same_type(const basic & other) const { return true; }
public:

	// substitutions
	ex subs(const exmap & m, unsigned options = 0) const { return inherited::subs(m, options); }

	// function mapping
	ex map(map_function & f) const { return inherited::map(f); }

	// degree/coeff
	int degree(const ex & s) const { return inherited::degree(s); }
	int ldegree(const ex & s) const { return inherited::ldegree(s); }
	ex coeff(const ex & s, int n = 1) const { return inherited::coeff(s, n); }

	// expand/collect
	ex expand(unsigned options = 0) const { return inherited::expand(options); }
	ex collect(const ex & s, bool distributed = false) const { return inherited::collect(s, distributed); }

	// differentiation and series expansion
protected:
	ex derivative(const symbol & s) const { return inherited::derivative(s); }
public:
	ex series(const relational & r, int order, unsigned options = 0) const { return inherited::series(r, order, options); }

	// rational functions
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0) const { return inherited::normal(repl, rev_lookup, level); }
	ex to_rational(exmap & repl) const { return inherited::to_rational(repl); }
	ex to_polynomial(exmap & repl) const { return inherited::to_polynomial(repl); }

	// polynomial algorithms
	numeric integer_content() const { return 1; }
	ex smod(const numeric & xi) const { return *this; }
	numeric max_coefficient() const { return 1; }

	// indexed objects
	exvector get_free_indices() const { return exvector(); }
	ex add_indexed(const ex & self, const ex & other) const { return self + other; }
	ex scalar_mul_indexed(const ex & self, const numeric & other) const { return self * ex(other); }
	bool contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const { return false; }

	// noncommutativity
	unsigned return_type() const { return return_types::commutative; }
	tinfo_t return_type_tinfo() const { return this; }

protected:
	bool is_equal_same_type(const basic & other) const
	{
		GINAC_ASSERT(is_a<structure>(other));
		const structure & o = static_cast<const structure &>(other);

		return struct_is_equal(&obj, &o.obj);
	}

	unsigned calchash() const { return inherited::calchash(); }

	// non-virtual functions in this class
public:
	// access to embedded structure
	const T *operator->() const { return &obj; }
	T &get_struct() { return obj; }
	const T &get_struct() const { return obj; }
private:
	T obj;
};


/** Default constructor */
template <class T, template <class> class CP>
structure<T, CP>::structure() : inherited(get_tinfo()) { }

/** Construct object from archive_node. */
template <class T, template <class> class CP>
structure<T, CP>::structure(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst) {}

/** Unarchive the object. */
template <class T, template <class> class CP>
ex structure<T, CP>::unarchive(const archive_node &n, lst &sym_lst)
{
	return (new structure(n, sym_lst))->setflag(status_flags::dynallocated);
}

/** Archive the object. */
template <class T, template <class> class CP>
void structure<T, CP>::archive(archive_node &n) const
{
	inherited::archive(n);
}

/** Compare two structures of the same type. */
template <class T, template <class> class CP>
int structure<T, CP>::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<structure>(other));
	const structure & o = static_cast<const structure &>(other);

	return struct_compare(&obj, &o.obj);
}

template <class T, template <class> class CP>
registered_class_info structure<T, CP>::reg_info = registered_class_info(registered_class_options(structure::get_class_name(), "basic", structure::next_structure_tinfo_key(), &structure::unarchive));


} // namespace GiNaC

#endif // ndef __GINAC_STRUCTURE_H__
