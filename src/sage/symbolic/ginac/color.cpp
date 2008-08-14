/** @file color.cpp
 *
 *  Implementation of GiNaC's color (SU(3) Lie algebra) objects. */

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

#include "color.h"
#include "idx.h"
#include "ncmul.h"
#include "symmetry.h"
#include "operators.h"
#include "numeric.h"
#include "mul.h"
#include "power.h" // for sqrt()
#include "symbol.h"
#include "archive.h"
#include "utils.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS(color, indexed)

const tinfo_static_t color::return_type_tinfo_static[256] = {{}};

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(su3one, tensor,
  print_func<print_dflt>(&su3one::do_print).
  print_func<print_latex>(&su3one::do_print_latex))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(su3t, tensor,
  print_func<print_dflt>(&su3t::do_print).
  print_func<print_latex>(&su3t::do_print))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(su3f, tensor,
  print_func<print_dflt>(&su3f::do_print).
  print_func<print_latex>(&su3f::do_print))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(su3d, tensor,
  print_func<print_dflt>(&su3d::do_print).
  print_func<print_latex>(&su3d::do_print))

//////////
// default constructors
//////////

color::color() : representation_label(0)
{
	tinfo_key = &color::tinfo_static;
}

DEFAULT_CTOR(su3one)
DEFAULT_CTOR(su3t)
DEFAULT_CTOR(su3f)
DEFAULT_CTOR(su3d)

//////////
// other constructors
//////////

/** Construct object without any color index. This constructor is for
 *  internal use only. Use the color_ONE() function instead.
 *  @see color_ONE */
color::color(const ex & b, unsigned char rl) : inherited(b), representation_label(rl)
{
	tinfo_key = &color::tinfo_static;
}

/** Construct object with one color index. This constructor is for internal
 *  use only. Use the color_T() function instead.
 *  @see color_T */
color::color(const ex & b, const ex & i1, unsigned char rl) : inherited(b, i1), representation_label(rl)
{
	tinfo_key = &color::tinfo_static;
}

color::color(unsigned char rl, const exvector & v, bool discardable) : inherited(not_symmetric(), v, discardable), representation_label(rl)
{
	tinfo_key = &color::tinfo_static;
}

color::color(unsigned char rl, std::auto_ptr<exvector> vp) : inherited(not_symmetric(), vp), representation_label(rl)
{
	tinfo_key = &color::tinfo_static;
}

//////////
// archiving
//////////

color::color(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	unsigned rl;
	n.find_unsigned("label", rl);
	representation_label = rl;
}

void color::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_unsigned("label", representation_label);
}

DEFAULT_UNARCHIVE(color)
DEFAULT_ARCHIVING(su3one)
DEFAULT_ARCHIVING(su3t)
DEFAULT_ARCHIVING(su3f)
DEFAULT_ARCHIVING(su3d)

//////////
// functions overriding virtual functions from base classes
//////////

int color::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<color>(other));
	const color &o = static_cast<const color &>(other);

	if (representation_label != o.representation_label) {
		// different representation label
		return representation_label < o.representation_label ? -1 : 1;
	}

	return inherited::compare_same_type(other);
}

bool color::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<color>(other));
	const color &o = static_cast<const color &>(other);

	return representation_label == o.representation_label;
}

DEFAULT_COMPARE(su3one)
DEFAULT_COMPARE(su3t)
DEFAULT_COMPARE(su3f)
DEFAULT_COMPARE(su3d)

DEFAULT_PRINT_LATEX(su3one, "ONE", "\\mathbb{1}")
DEFAULT_PRINT(su3t, "T")
DEFAULT_PRINT(su3f, "f")
DEFAULT_PRINT(su3d, "d")

/** Perform automatic simplification on noncommutative product of color
 *  objects. This removes superfluous ONEs. */
ex color::eval_ncmul(const exvector & v) const
{
	exvector s;
	s.reserve(v.size());

	// Remove superfluous ONEs
	exvector::const_iterator it = v.begin(), itend = v.end();
	while (it != itend) {
		if (!is_a<su3one>(it->op(0)))
			s.push_back(*it);
		it++;
	}

	if (s.empty())
		return color(su3one(), representation_label);
	else
		return hold_ncmul(s);
}

ex color::thiscontainer(const exvector & v) const
{
	return color(representation_label, v);
}

ex color::thiscontainer(std::auto_ptr<exvector> vp) const
{
	return color(representation_label, vp);
}

/** Given a vector iv3 of three indices and a vector iv2 of two indices that
 *  is a subset of iv3, return the (free) index that is in iv3 but not in
 *  iv2 and the sign introduced by permuting that index to the front.
 *
 *  @param iv3 Vector of 3 indices
 *  @param iv2 Vector of 2 indices, must be a subset of iv3
 *  @param sig Returs sign introduced by index permutation
 *  @return the free index (the one that is in iv3 but not in iv2) */
static ex permute_free_index_to_front(const exvector & iv3, const exvector & iv2, int & sig)
{
	GINAC_ASSERT(iv3.size() == 3);
	GINAC_ASSERT(iv2.size() == 2);

	sig = 1;

#define TEST_PERMUTATION(A,B,C,P) \
	if (iv3[B].is_equal(iv2[0]) && iv3[C].is_equal(iv2[1])) { \
		sig = P; \
		return iv3[A]; \
	}
	
	TEST_PERMUTATION(0,1,2,  1);
	TEST_PERMUTATION(0,2,1, -1);
	TEST_PERMUTATION(1,0,2, -1);
	TEST_PERMUTATION(1,2,0,  1);
	TEST_PERMUTATION(2,0,1,  1);
	TEST_PERMUTATION(2,1,0, -1);

	throw(std::logic_error("permute_free_index_to_front(): no valid permutation found"));
}

/** Automatic symbolic evaluation of indexed symmetric structure constant. */
ex su3d::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(i.nops() == 4);
	GINAC_ASSERT(is_a<su3d>(i.op(0)));

	// Convolutions are zero
	if (!(static_cast<const indexed &>(i).get_dummy_indices().empty()))
		return _ex0;

	// Numeric evaluation
	if (static_cast<const indexed &>(i).all_index_values_are(info_flags::nonnegint)) {

		// Sort indices
		int v[3];
		for (unsigned j=0; j<3; j++)
			v[j] = ex_to<numeric>(ex_to<idx>(i.op(j + 1)).get_value()).to_int();
		if (v[0] > v[1]) std::swap(v[0], v[1]);
		if (v[0] > v[2]) std::swap(v[0], v[2]);
		if (v[1] > v[2]) std::swap(v[1], v[2]);

#define CMPINDICES(A,B,C) ((v[0] == (A)) && (v[1] == (B)) && (v[2] == (C)))

		// Check for non-zero elements
		if (CMPINDICES(1,4,6) || CMPINDICES(1,5,7) || CMPINDICES(2,5,6)
		 || CMPINDICES(3,4,4) || CMPINDICES(3,5,5))
			return _ex1_2;
		else if (CMPINDICES(2,4,7) || CMPINDICES(3,6,6) || CMPINDICES(3,7,7))
			return _ex_1_2;
		else if (CMPINDICES(1,1,8) || CMPINDICES(2,2,8) || CMPINDICES(3,3,8))
			return sqrt(_ex3)*_ex1_3;
		else if (CMPINDICES(8,8,8))
			return sqrt(_ex3)*_ex_1_3;
		else if (CMPINDICES(4,4,8) || CMPINDICES(5,5,8)
		      || CMPINDICES(6,6,8) || CMPINDICES(7,7,8))
			return sqrt(_ex3)/_ex_6;
		else
			return _ex0;
	}

	// No further simplifications
	return i.hold();
}

/** Automatic symbolic evaluation of indexed antisymmetric structure constant. */
ex su3f::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(i.nops() == 4);
	GINAC_ASSERT(is_a<su3f>(i.op(0)));

	// Numeric evaluation
	if (static_cast<const indexed &>(i).all_index_values_are(info_flags::nonnegint)) {

		// Sort indices, remember permutation sign
		int v[3];
		for (unsigned j=0; j<3; j++)
			v[j] = ex_to<numeric>(ex_to<idx>(i.op(j + 1)).get_value()).to_int();
		int sign = 1;
		if (v[0] > v[1]) { std::swap(v[0], v[1]); sign = -sign; }
		if (v[0] > v[2]) { std::swap(v[0], v[2]); sign = -sign; }
		if (v[1] > v[2]) { std::swap(v[1], v[2]); sign = -sign; }

		// Check for non-zero elements
		if (CMPINDICES(1,2,3))
			return sign;
		else if (CMPINDICES(1,4,7) || CMPINDICES(2,4,6)
		      || CMPINDICES(2,5,7) || CMPINDICES(3,4,5))
			return _ex1_2 * sign;
		else if (CMPINDICES(1,5,6) || CMPINDICES(3,6,7))
			return _ex_1_2 * sign;
		else if (CMPINDICES(4,5,8) || CMPINDICES(6,7,8))
			return sqrt(_ex3)/2 * sign;
		else
			return _ex0;
	}

	// No further simplifications
	return i.hold();
}


/** Contraction of generator with something else. */
bool su3t::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(self->nops() == 2);
	GINAC_ASSERT(is_a<su3t>(self->op(0)));
	unsigned char rl = ex_to<color>(*self).get_representation_label();

	if (is_exactly_a<su3t>(other->op(0))) {

		// Contraction only makes sense if the represenation labels are equal
		GINAC_ASSERT(is_a<color>(*other));
		if (ex_to<color>(*other).get_representation_label() != rl)
			return false;

		// T.a T.a = 4/3 ONE
		if (other - self == 1) {
			*self = numeric(4, 3);
			*other = color_ONE(rl);
			return true;

		// T.a T.b T.a = -1/6 T.b
		} else if (other - self == 2
		        && is_a<color>(self[1])) {
			*self = numeric(-1, 6);
			*other = _ex1;
			return true;

		// T.a S T.a = 1/2 Tr(S) - 1/6 S
		} else {
			exvector::iterator it = self + 1;
			while (it != other) {
				if (!is_a<color>(*it)) {
					return false;
				}
				it++;
			}

			it = self + 1;
			ex S = _ex1;
			while (it != other) {
				S *= *it;
				*it++ = _ex1;
			}

			*self = color_trace(S, rl) * color_ONE(rl) / 2 - S / 6;
			*other = _ex1;
			return true;
		}
	}

	return false;
}

/** Contraction of an indexed symmetric structure constant with something else. */
bool su3d::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(self->nops() == 4);
	GINAC_ASSERT(is_a<su3d>(self->op(0)));

	if (is_exactly_a<su3d>(other->op(0))) {

		// Find the dummy indices of the contraction
		exvector self_indices = ex_to<indexed>(*self).get_indices();
		exvector other_indices = ex_to<indexed>(*other).get_indices();
		exvector all_indices = self_indices;
		all_indices.insert(all_indices.end(), other_indices.begin(), other_indices.end());
		exvector free_indices, dummy_indices;
		find_free_and_dummy(all_indices, free_indices, dummy_indices);

		// d.abc d.abc = 40/3
		if (dummy_indices.size() == 3) {
			*self = numeric(40, 3);
			*other = _ex1;
			return true;

		// d.akl d.bkl = 5/3 delta.ab
		} else if (dummy_indices.size() == 2) {
			exvector a;
			std::back_insert_iterator<exvector> ita(a);
			ita = set_difference(self_indices.begin(), self_indices.end(), dummy_indices.begin(), dummy_indices.end(), ita, ex_is_less());
			ita = set_difference(other_indices.begin(), other_indices.end(), dummy_indices.begin(), dummy_indices.end(), ita, ex_is_less());
			GINAC_ASSERT(a.size() == 2);
			*self = numeric(5, 3) * delta_tensor(a[0], a[1]);
			*other = _ex1;
			return true;
		}

	} else if (is_exactly_a<su3t>(other->op(0))) {

		// d.abc T.b T.c = 5/6 T.a
		if (other+1 != v.end()
		 && is_exactly_a<su3t>(other[1].op(0))
		 && ex_to<indexed>(*self).has_dummy_index_for(other[1].op(1))) {

			exvector self_indices = ex_to<indexed>(*self).get_indices();
			exvector dummy_indices;
			dummy_indices.push_back(other[0].op(1));
			dummy_indices.push_back(other[1].op(1));
			int sig;
			ex a = permute_free_index_to_front(self_indices, dummy_indices, sig);
			*self = numeric(5, 6);
			other[0] = color_T(a, ex_to<color>(other[0]).get_representation_label());
			other[1] = _ex1;
			return true;
		}
	}

	return false;
}

/** Contraction of an indexed antisymmetric structure constant with something else. */
bool su3f::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(self->nops() == 4);
	GINAC_ASSERT(is_a<su3f>(self->op(0)));

	if (is_exactly_a<su3f>(other->op(0))) { // f*d is handled by su3d class

		// Find the dummy indices of the contraction
		exvector dummy_indices;
		dummy_indices = ex_to<indexed>(*self).get_dummy_indices(ex_to<indexed>(*other));

		// f.abc f.abc = 24
		if (dummy_indices.size() == 3) {
			*self = 24;
			*other = _ex1;
			return true;

		// f.akl f.bkl = 3 delta.ab
		} else if (dummy_indices.size() == 2) {
			int sign1, sign2;
			ex a = permute_free_index_to_front(ex_to<indexed>(*self).get_indices(), dummy_indices, sign1);
			ex b = permute_free_index_to_front(ex_to<indexed>(*other).get_indices(), dummy_indices, sign2);
			*self = sign1 * sign2 * 3 * delta_tensor(a, b);
			*other = _ex1;
			return true;
		}

	} else if (is_exactly_a<su3t>(other->op(0))) {

		// f.abc T.b T.c = 3/2 I T.a
		if (other+1 != v.end()
		 && is_exactly_a<su3t>(other[1].op(0))
		 && ex_to<indexed>(*self).has_dummy_index_for(other[1].op(1))) {

			exvector self_indices = ex_to<indexed>(*self).get_indices();
			exvector dummy_indices;
			dummy_indices.push_back(other[0].op(1));
			dummy_indices.push_back(other[1].op(1));
			int sig;
			ex a = permute_free_index_to_front(self_indices, dummy_indices, sig);
			*self = numeric(3, 2) * sig * I;
			other[0] = color_T(a, ex_to<color>(other[0]).get_representation_label());
			other[1] = _ex1;
			return true;
		}
	}

	return false;
}

//////////
// global functions
//////////

ex color_ONE(unsigned char rl)
{
	static ex ONE = (new su3one)->setflag(status_flags::dynallocated);
	return color(ONE, rl);
}

ex color_T(const ex & a, unsigned char rl)
{
	static ex t = (new su3t)->setflag(status_flags::dynallocated);

	if (!is_a<idx>(a))
		throw(std::invalid_argument("indices of color_T must be of type idx"));
	if (!ex_to<idx>(a).get_dim().is_equal(8))
		throw(std::invalid_argument("index dimension for color_T must be 8"));

	return color(t, a, rl);
}

ex color_f(const ex & a, const ex & b, const ex & c)
{
	static ex f = (new su3f)->setflag(status_flags::dynallocated);

	if (!is_a<idx>(a) || !is_a<idx>(b) || !is_a<idx>(c))
		throw(std::invalid_argument("indices of color_f must be of type idx"));
	if (!ex_to<idx>(a).get_dim().is_equal(8) || !ex_to<idx>(b).get_dim().is_equal(8) || !ex_to<idx>(c).get_dim().is_equal(8))
		throw(std::invalid_argument("index dimension for color_f must be 8"));

	return indexed(f, antisymmetric3(), a, b, c);
}

ex color_d(const ex & a, const ex & b, const ex & c)
{
	static ex d = (new su3d)->setflag(status_flags::dynallocated);

	if (!is_a<idx>(a) || !is_a<idx>(b) || !is_a<idx>(c))
		throw(std::invalid_argument("indices of color_d must be of type idx"));
	if (!ex_to<idx>(a).get_dim().is_equal(8) || !ex_to<idx>(b).get_dim().is_equal(8) || !ex_to<idx>(c).get_dim().is_equal(8))
		throw(std::invalid_argument("index dimension for color_d must be 8"));

	return indexed(d, symmetric3(), a, b, c);
}

ex color_h(const ex & a, const ex & b, const ex & c)
{
	return color_d(a, b, c) + I * color_f(a, b, c);
}

/** Check whether a given tinfo key (as returned by return_type_tinfo()
 *  is that of a color object (with an arbitrary representation label). */
static bool is_color_tinfo(tinfo_t ti)
{
	p_int start_loc=(p_int)&color::return_type_tinfo_static;
	return (p_int)ti>=start_loc && (p_int)ti<start_loc+256;
}

/** Extract representation label from tinfo key (as returned by
 *  return_type_tinfo()). */
static unsigned char get_representation_label(tinfo_t ti)
{
	return (unsigned char)((p_int)ti-(p_int)&color::return_type_tinfo_static);
}

ex color_trace(const ex & e, const std::set<unsigned char> & rls)
{
	if (is_a<color>(e)) {

		unsigned char rl = ex_to<color>(e).get_representation_label();

		// Are we taking the trace over this object's representation label?
		if (rls.find(rl) == rls.end())
			return e;

		// Yes, all generators are traceless, except for color_ONE
		if (is_a<su3one>(e.op(0)))
			return _ex3;
		else
			return _ex0;

	} else if (is_exactly_a<mul>(e)) {

		// Trace of product: pull out non-color factors
		ex prod = _ex1;
		for (size_t i=0; i<e.nops(); i++) {
			const ex &o = e.op(i);
			if (is_color_tinfo(o.return_type_tinfo()))
				prod *= color_trace(o, rls);
			else
				prod *= o;
		}
		return prod;

	} else if (is_exactly_a<ncmul>(e)) {

		unsigned char rl = get_representation_label(e.return_type_tinfo());

		// Are we taking the trace over this string's representation label?
		if (rls.find(rl) == rls.end())
			return e;

		// Yes, expand product if necessary
		ex e_expanded = e.expand();
		if (!is_a<ncmul>(e_expanded))
			return color_trace(e_expanded, rls);

		size_t num = e.nops();

		if (num == 2) {

			// Tr T_a T_b = 1/2 delta_a_b
			return delta_tensor(e.op(0).op(1), e.op(1).op(1)) / 2;

		} else if (num == 3) {

			// Tr T_a T_b T_c = 1/4 h_a_b_c
			return color_h(e.op(0).op(1), e.op(1).op(1), e.op(2).op(1)) / 4;

		} else {

			// Traces of 4 or more generators are computed recursively:
			// Tr T_a1 .. T_an =
			//     1/6 delta_a(n-1)_an Tr T_a1 .. T_a(n-2)
			//   + 1/2 h_a(n-1)_an_k Tr T_a1 .. T_a(n-2) T_k
			const ex &last_index = e.op(num - 1).op(1);
			const ex &next_to_last_index = e.op(num - 2).op(1);
			idx summation_index((new symbol)->setflag(status_flags::dynallocated), 8);

			exvector v1;
			v1.reserve(num - 2);
			for (size_t i=0; i<num-2; i++)
				v1.push_back(e.op(i));

			exvector v2 = v1;
			v2.push_back(color_T(summation_index, rl));

			return delta_tensor(next_to_last_index, last_index) * color_trace(ncmul(v1), rl) / 6
			       + color_h(next_to_last_index, last_index, summation_index) * color_trace(ncmul(v2), rl) / 2;
		}

	} else if (e.nops() > 0) {

		// Trace maps to all other container classes (this includes sums)
		pointer_to_map_function_1arg<const std::set<unsigned char> &> fcn(color_trace, rls);
		return e.map(fcn);

	} else
		return _ex0;
}

ex color_trace(const ex & e, const lst & rll)
{
	// Convert list to set
	std::set<unsigned char> rls;
	for (lst::const_iterator i = rll.begin(); i != rll.end(); ++i) {
		if (i->info(info_flags::nonnegint))
			rls.insert(ex_to<numeric>(*i).to_int());
	}

	return color_trace(e, rls);
}

ex color_trace(const ex & e, unsigned char rl)
{
	// Convert label to set
	std::set<unsigned char> rls;
	rls.insert(rl);

	return color_trace(e, rls);
}

} // namespace GiNaC
