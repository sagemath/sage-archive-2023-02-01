/** @file tensor.cpp
 *
 *  Implementation of GiNaC's special tensors. */

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
#include <vector>

#include "tensor.h"
#include "idx.h"
#include "indexed.h"
#include "symmetry.h"
#include "relational.h"
#include "operators.h"
#include "lst.h"
#include "numeric.h"
#include "matrix.h"
#include "archive.h"
#include "utils.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS(tensor, basic)

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(tensdelta, tensor,
  print_func<print_dflt>(&tensdelta::do_print).
  print_func<print_latex>(&tensdelta::do_print_latex))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(tensmetric, tensor,
  print_func<print_dflt>(&tensmetric::do_print).
  print_func<print_latex>(&tensmetric::do_print))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(minkmetric, tensmetric,
  print_func<print_dflt>(&minkmetric::do_print).
  print_func<print_latex>(&minkmetric::do_print_latex))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(spinmetric, tensmetric,
  print_func<print_dflt>(&spinmetric::do_print).
  print_func<print_latex>(&spinmetric::do_print_latex))

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(tensepsilon, tensor,
  print_func<print_dflt>(&tensepsilon::do_print).
  print_func<print_latex>(&tensepsilon::do_print_latex))

//////////
// constructors
//////////

tensor::tensor() : inherited(&tensor::tinfo_static)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

DEFAULT_CTOR(tensdelta)
DEFAULT_CTOR(tensmetric)

minkmetric::minkmetric() : pos_sig(false)
{
	tinfo_key = &minkmetric::tinfo_static;
}

spinmetric::spinmetric()
{
	tinfo_key = &spinmetric::tinfo_static;
}

minkmetric::minkmetric(bool ps) : pos_sig(ps)
{
	tinfo_key = &minkmetric::tinfo_static;
}

tensepsilon::tensepsilon() : minkowski(false), pos_sig(false)
{
	tinfo_key = &tensepsilon::tinfo_static;
}

tensepsilon::tensepsilon(bool mink, bool ps) : minkowski(mink), pos_sig(ps)
{
	tinfo_key = &tensepsilon::tinfo_static;
}

//////////
// archiving
//////////

DEFAULT_ARCHIVING(tensor)
DEFAULT_ARCHIVING(tensdelta)
DEFAULT_ARCHIVING(tensmetric)
DEFAULT_ARCHIVING(spinmetric)
DEFAULT_UNARCHIVE(minkmetric)
DEFAULT_UNARCHIVE(tensepsilon)

minkmetric::minkmetric(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	n.find_bool("pos_sig", pos_sig);
}

void minkmetric::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_bool("pos_sig", pos_sig);
}

tensepsilon::tensepsilon(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	n.find_bool("minkowski", minkowski);
	n.find_bool("pos_sig", pos_sig);
}

void tensepsilon::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_bool("minkowski", minkowski);
	n.add_bool("pos_sig", pos_sig);
}

//////////
// functions overriding virtual functions from base classes
//////////

DEFAULT_COMPARE(tensor)
DEFAULT_COMPARE(tensdelta)
DEFAULT_COMPARE(tensmetric)
DEFAULT_COMPARE(spinmetric)

bool tensdelta::info(unsigned inf) const
{
	if(inf == info_flags::real)
		return true;

	return false;
}

bool tensmetric::info(unsigned inf) const
{
	if(inf == info_flags::real)
		return true;

	return false;
}

int minkmetric::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<minkmetric>(other));
	const minkmetric &o = static_cast<const minkmetric &>(other);

	if (pos_sig != o.pos_sig)
		return pos_sig ? -1 : 1;
	else
		return inherited::compare_same_type(other);
}

bool minkmetric::info(unsigned inf) const
{
	if(inf == info_flags::real)
		return true;

	return false;
}

int tensepsilon::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<tensepsilon>(other));
	const tensepsilon &o = static_cast<const tensepsilon &>(other);

	if (minkowski != o.minkowski)
		return minkowski ? -1 : 1;
	else if (pos_sig != o.pos_sig)
		return pos_sig ? -1 : 1;
	else
		return inherited::compare_same_type(other);
}

bool tensepsilon::info(unsigned inf) const
{
	if(inf == info_flags::real)
		return true;

	return false;
}

bool spinmetric::info(unsigned inf) const
{
	if(inf == info_flags::real)
		return true;

	return false;
}

DEFAULT_PRINT_LATEX(tensdelta, "delta", "\\delta")
DEFAULT_PRINT(tensmetric, "g")
DEFAULT_PRINT_LATEX(minkmetric, "eta", "\\eta")
DEFAULT_PRINT_LATEX(spinmetric, "eps", "\\varepsilon")
DEFAULT_PRINT_LATEX(tensepsilon, "eps", "\\varepsilon")

/** Automatic symbolic evaluation of an indexed delta tensor. */
ex tensdelta::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(i.nops() == 3);
	GINAC_ASSERT(is_a<tensdelta>(i.op(0)));

	const idx & i1 = ex_to<idx>(i.op(1));
	const idx & i2 = ex_to<idx>(i.op(2));

	// The dimension of the indices must be equal, otherwise we use the minimal
	// dimension
	if (!i1.get_dim().is_equal(i2.get_dim())) {
		ex min_dim = i1.minimal_dim(i2);
		exmap m;
		m[i1] = i1.replace_dim(min_dim);
		m[i2] = i2.replace_dim(min_dim);
		return i.subs(m, subs_options::no_pattern);
	}

	// Trace of delta tensor is the (effective) dimension of the space
	if (is_dummy_pair(i1, i2)) {
		try {
			return i1.minimal_dim(i2);
		} catch (std::exception &e) {
			return i.hold();
		}
	}

	// Numeric evaluation
	if (static_cast<const indexed &>(i).all_index_values_are(info_flags::integer)) {
		int n1 = ex_to<numeric>(i1.get_value()).to_int(), n2 = ex_to<numeric>(i2.get_value()).to_int();
		if (n1 == n2)
			return _ex1;
		else
			return _ex0;
	}

	// No further simplifications
	return i.hold();
}

/** Automatic symbolic evaluation of an indexed metric tensor. */
ex tensmetric::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(i.nops() == 3);
	GINAC_ASSERT(is_a<tensmetric>(i.op(0)));
	GINAC_ASSERT(is_a<varidx>(i.op(1)));
	GINAC_ASSERT(is_a<varidx>(i.op(2)));

	const varidx & i1 = ex_to<varidx>(i.op(1));
	const varidx & i2 = ex_to<varidx>(i.op(2));

	// The dimension of the indices must be equal, otherwise we use the minimal
	// dimension
	if (!i1.get_dim().is_equal(i2.get_dim())) {
		ex min_dim = i1.minimal_dim(i2);
		exmap m;
		m[i1] = i1.replace_dim(min_dim);
		m[i2] = i2.replace_dim(min_dim);
		return i.subs(m, subs_options::no_pattern);
	}

	// A metric tensor with one covariant and one contravariant index gets
	// replaced by a delta tensor
	if (i1.is_covariant() != i2.is_covariant())
		return delta_tensor(i1, i2);

	// No further simplifications
	return i.hold();
}

/** Automatic symbolic evaluation of an indexed Lorentz metric tensor. */
ex minkmetric::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(i.nops() == 3);
	GINAC_ASSERT(is_a<minkmetric>(i.op(0)));
	GINAC_ASSERT(is_a<varidx>(i.op(1)));
	GINAC_ASSERT(is_a<varidx>(i.op(2)));

	const varidx & i1 = ex_to<varidx>(i.op(1));
	const varidx & i2 = ex_to<varidx>(i.op(2));

	// Numeric evaluation
	if (static_cast<const indexed &>(i).all_index_values_are(info_flags::nonnegint)) {
		int n1 = ex_to<numeric>(i1.get_value()).to_int(), n2 = ex_to<numeric>(i2.get_value()).to_int();
		if (n1 != n2)
			return _ex0;
		else if (n1 == 0)
			return pos_sig ? _ex_1 : _ex1;
		else
			return pos_sig ? _ex1 : _ex_1;
	}

	// Perform the usual evaluations of a metric tensor
	return inherited::eval_indexed(i);
}

/** Automatic symbolic evaluation of an indexed metric tensor. */
ex spinmetric::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(i.nops() == 3);
	GINAC_ASSERT(is_a<spinmetric>(i.op(0)));
	GINAC_ASSERT(is_a<spinidx>(i.op(1)));
	GINAC_ASSERT(is_a<spinidx>(i.op(2)));

	const spinidx & i1 = ex_to<spinidx>(i.op(1));
	const spinidx & i2 = ex_to<spinidx>(i.op(2));

	// Convolutions are zero
	if (!(static_cast<const indexed &>(i).get_dummy_indices().empty()))
		return _ex0;

	// Numeric evaluation
	if (static_cast<const indexed &>(i).all_index_values_are(info_flags::nonnegint)) {
		int n1 = ex_to<numeric>(i1.get_value()).to_int(), n2 = ex_to<numeric>(i2.get_value()).to_int();
		if (n1 == n2)
			return _ex0;
		else if (n1 < n2)
			return _ex1;
		else
			return _ex_1;
	}

	// No further simplifications
	return i.hold();
}

/** Automatic symbolic evaluation of an indexed epsilon tensor. */
ex tensepsilon::eval_indexed(const basic & i) const
{
	GINAC_ASSERT(is_a<indexed>(i));
	GINAC_ASSERT(i.nops() > 1);
	GINAC_ASSERT(is_a<tensepsilon>(i.op(0)));

	// Convolutions are zero
	if (!(static_cast<const indexed &>(i).get_dummy_indices().empty()))
		return _ex0;

	// Numeric evaluation
	if (static_cast<const indexed &>(i).all_index_values_are(info_flags::nonnegint)) {

		// Get sign of index permutation (the indices should already be in
		// a canonic order but we can't assume what exactly that order is)
		std::vector<int> v;
		v.reserve(i.nops() - 1);
		for (size_t j=1; j<i.nops(); j++)
			v.push_back(ex_to<numeric>(ex_to<idx>(i.op(j)).get_value()).to_int());
		int sign = permutation_sign(v.begin(), v.end());

		// In a Minkowski space, check for covariant indices
		if (minkowski) {
			for (size_t j=1; j<i.nops(); j++) {
				const ex & x = i.op(j);
				if (!is_a<varidx>(x))
					throw(std::runtime_error("indices of epsilon tensor in Minkowski space must be of type varidx"));
				if (ex_to<varidx>(x).is_covariant())
					if (ex_to<idx>(x).get_value().is_zero())
						sign = (pos_sig ? -sign : sign);
					else
						sign = (pos_sig ? sign : -sign);
			}
		}

		return sign;
	}

	// No further simplifications
	return i.hold();
}

bool tensor::replace_contr_index(exvector::iterator self, exvector::iterator other) const
{
	// Try to contract the first index
	const idx *self_idx = &ex_to<idx>(self->op(1));
	const idx *free_idx = &ex_to<idx>(self->op(2));
	bool first_index_tried = false;

again:
	if (self_idx->is_symbolic()) {
		for (size_t i=1; i<other->nops(); i++) {
			if (! is_a<idx>(other->op(i)))
				continue;
			const idx &other_idx = ex_to<idx>(other->op(i));
			if (is_dummy_pair(*self_idx, other_idx)) {

				// Contraction found, remove this tensor and substitute the
				// index in the second object
				try {
					// minimal_dim() throws an exception when index dimensions are not comparable
					ex min_dim = self_idx->minimal_dim(other_idx);
					*other = other->subs(other_idx == free_idx->replace_dim(min_dim));
					*self = _ex1; // *other is assigned first because assigning *self invalidates free_idx
					return true;
				} catch (std::exception &e) {
					return false;
				}
			}
		}
	}

	if (!first_index_tried) {

		// No contraction with the first index found, try the second index
		self_idx = &ex_to<idx>(self->op(2));
		free_idx = &ex_to<idx>(self->op(1));
		first_index_tried = true;
		goto again;
	}

	return false;
}

/** Contraction of an indexed delta tensor with something else. */
bool tensdelta::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(self->nops() == 3);
	GINAC_ASSERT(is_a<tensdelta>(self->op(0)));

	// Replace the dummy index with this tensor's other index and remove
	// the tensor (this is valid for contractions with all other tensors)
	return replace_contr_index(self, other);
}

/** Contraction of an indexed metric tensor with something else. */
bool tensmetric::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(self->nops() == 3);
	GINAC_ASSERT(is_a<tensmetric>(self->op(0)));

	// If contracting with the delta tensor, let the delta do it
	// (don't raise/lower delta indices)
	if (is_a<tensdelta>(other->op(0)))
		return false;

	// Replace the dummy index with this tensor's other index and remove
	// the tensor
	return replace_contr_index(self, other);
}

/** Contraction of an indexed spinor metric with something else. */
bool spinmetric::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(self->nops() == 3);
	GINAC_ASSERT(is_a<spinmetric>(self->op(0)));

	// Contractions between spinor metrics
	if (is_a<spinmetric>(other->op(0))) {
		const idx &self_i1 = ex_to<idx>(self->op(1));
		const idx &self_i2 = ex_to<idx>(self->op(2));
		const idx &other_i1 = ex_to<idx>(other->op(1));
		const idx &other_i2 = ex_to<idx>(other->op(2));

		if (is_dummy_pair(self_i1, other_i1)) {
			if (is_dummy_pair(self_i2, other_i2))
				*self = _ex2;
			else
				*self = delta_tensor(self_i2, other_i2);
			*other = _ex1;
			return true;
		} else if (is_dummy_pair(self_i1, other_i2)) {
			if (is_dummy_pair(self_i2, other_i1))
				*self = _ex_2;
			else
				*self = -delta_tensor(self_i2, other_i1);
			*other = _ex1;
			return true;
		} else if (is_dummy_pair(self_i2, other_i1)) {
			*self = -delta_tensor(self_i1, other_i2);
			*other = _ex1;
			return true;
		} else if (is_dummy_pair(self_i2, other_i2)) {
			*self = delta_tensor(self_i1, other_i1);
			*other = _ex1;
			return true;
		}
	}

	// If contracting with the delta tensor, let the delta do it
	// (don't raise/lower delta indices)
	if (is_a<tensdelta>(other->op(0)))
		return false;

	// Try to contract first index
	const idx *self_idx = &ex_to<idx>(self->op(1));
	const idx *free_idx = &ex_to<idx>(self->op(2));
	bool first_index_tried = false;
	int sign = 1;

again:
	if (self_idx->is_symbolic()) {
		for (size_t i=1; i<other->nops(); i++) {
			const idx &other_idx = ex_to<idx>(other->op(i));
			if (is_dummy_pair(*self_idx, other_idx)) {

				// Contraction found, remove metric tensor and substitute
				// index in second object (assign *self last because this
				// invalidates free_idx)
				*other = other->subs(other_idx == *free_idx);
				*self = (static_cast<const spinidx *>(self_idx)->is_covariant() ? sign : -sign);
				return true;
			}
		}
	}

	if (!first_index_tried) {

		// No contraction with first index found, try second index
		self_idx = &ex_to<idx>(self->op(2));
		free_idx = &ex_to<idx>(self->op(1));
		first_index_tried = true;
		sign = -sign;
		goto again;
	}

	return false;
}

/** Contraction of epsilon tensor with something else. */
bool tensepsilon::contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const
{
	GINAC_ASSERT(is_a<indexed>(*self));
	GINAC_ASSERT(is_a<indexed>(*other));
	GINAC_ASSERT(is_a<tensepsilon>(self->op(0)));
	size_t num = self->nops() - 1;

	if (is_exactly_a<tensepsilon>(other->op(0)) && num+1 == other->nops()) {

		// Contraction of two epsilon tensors is a determinant
		bool variance = is_a<varidx>(self->op(1));
		matrix M(num, num);
		for (size_t i=0; i<num; i++) {
			for (size_t j=0; j<num; j++) {
				if (minkowski)
					M(i, j) = lorentz_g(self->op(i+1), other->op(j+1), pos_sig);
				else if (variance)
					M(i, j) = metric_tensor(self->op(i+1), other->op(j+1));
				else
					M(i, j) = delta_tensor(self->op(i+1), other->op(j+1));
			}
		}
		int sign = minkowski ? -1 : 1;
		*self = sign * M.determinant().simplify_indexed();
		*other = _ex1;
		return true;
	}

	return false;
}

//////////
// global functions
//////////

ex delta_tensor(const ex & i1, const ex & i2)
{
	static ex delta = (new tensdelta)->setflag(status_flags::dynallocated);

	if (!is_a<idx>(i1) || !is_a<idx>(i2))
		throw(std::invalid_argument("indices of delta tensor must be of type idx"));

	return indexed(delta, symmetric2(), i1, i2);
}

ex metric_tensor(const ex & i1, const ex & i2)
{
	static ex metric = (new tensmetric)->setflag(status_flags::dynallocated);

	if (!is_a<varidx>(i1) || !is_a<varidx>(i2))
		throw(std::invalid_argument("indices of metric tensor must be of type varidx"));

	return indexed(metric, symmetric2(), i1, i2);
}

ex lorentz_g(const ex & i1, const ex & i2, bool pos_sig)
{
	static ex metric_neg = (new minkmetric(false))->setflag(status_flags::dynallocated);
	static ex metric_pos = (new minkmetric(true))->setflag(status_flags::dynallocated);

	if (!is_a<varidx>(i1) || !is_a<varidx>(i2))
		throw(std::invalid_argument("indices of metric tensor must be of type varidx"));

	return indexed(pos_sig ? metric_pos : metric_neg, symmetric2(), i1, i2);
}

ex spinor_metric(const ex & i1, const ex & i2)
{
	static ex metric = (new spinmetric)->setflag(status_flags::dynallocated);

	if (!is_a<spinidx>(i1) || !is_a<spinidx>(i2))
		throw(std::invalid_argument("indices of spinor metric must be of type spinidx"));
	if (!ex_to<idx>(i1).get_dim().is_equal(2) || !ex_to<idx>(i2).get_dim().is_equal(2))
		throw(std::runtime_error("index dimension for spinor metric must be 2"));

	return indexed(metric, antisymmetric2(), i1, i2);
}

ex epsilon_tensor(const ex & i1, const ex & i2)
{
	static ex epsilon = (new tensepsilon)->setflag(status_flags::dynallocated);

	if (!is_a<idx>(i1) || !is_a<idx>(i2))
		throw(std::invalid_argument("indices of epsilon tensor must be of type idx"));

	ex dim = ex_to<idx>(i1).get_dim();
	if (!dim.is_equal(ex_to<idx>(i2).get_dim()))
		throw(std::invalid_argument("all indices of epsilon tensor must have the same dimension"));
	if (!ex_to<idx>(i1).get_dim().is_equal(_ex2))
		throw(std::runtime_error("index dimension of epsilon tensor must match number of indices"));

	if(is_a<wildcard>(i1.op(0))||is_a<wildcard>(i2.op(0)))
		return indexed(epsilon, antisymmetric2(), i1, i2).hold();

	return indexed(epsilon, antisymmetric2(), i1, i2);
}

ex epsilon_tensor(const ex & i1, const ex & i2, const ex & i3)
{
	static ex epsilon = (new tensepsilon)->setflag(status_flags::dynallocated);

	if (!is_a<idx>(i1) || !is_a<idx>(i2) || !is_a<idx>(i3))
		throw(std::invalid_argument("indices of epsilon tensor must be of type idx"));

	ex dim = ex_to<idx>(i1).get_dim();
	if (!dim.is_equal(ex_to<idx>(i2).get_dim()) || !dim.is_equal(ex_to<idx>(i3).get_dim()))
		throw(std::invalid_argument("all indices of epsilon tensor must have the same dimension"));
	if (!ex_to<idx>(i1).get_dim().is_equal(_ex3))
		throw(std::runtime_error("index dimension of epsilon tensor must match number of indices"));

	if(is_a<wildcard>(i1.op(0))||is_a<wildcard>(i2.op(0))||is_a<wildcard>(i3.op(0)))
		return indexed(epsilon, antisymmetric3(), i1, i2, i3).hold();

	return indexed(epsilon, antisymmetric3(), i1, i2, i3);
}

ex lorentz_eps(const ex & i1, const ex & i2, const ex & i3, const ex & i4, bool pos_sig)
{
	static ex epsilon_neg = (new tensepsilon(true, false))->setflag(status_flags::dynallocated);
	static ex epsilon_pos = (new tensepsilon(true, true))->setflag(status_flags::dynallocated);

	if (!is_a<varidx>(i1) || !is_a<varidx>(i2) || !is_a<varidx>(i3) || !is_a<varidx>(i4))
		throw(std::invalid_argument("indices of Lorentz epsilon tensor must be of type varidx"));

	ex dim = ex_to<idx>(i1).get_dim();
	if (!dim.is_equal(ex_to<idx>(i2).get_dim()) || !dim.is_equal(ex_to<idx>(i3).get_dim()) || !dim.is_equal(ex_to<idx>(i4).get_dim()))
		throw(std::invalid_argument("all indices of epsilon tensor must have the same dimension"));
	if (!ex_to<idx>(i1).get_dim().is_equal(_ex4))
		throw(std::runtime_error("index dimension of epsilon tensor must match number of indices"));

	if(is_a<wildcard>(i1.op(0))||is_a<wildcard>(i2.op(0))||is_a<wildcard>(i3.op(0))||is_a<wildcard>(i4.op(0)))
		return indexed(pos_sig ? epsilon_pos : epsilon_neg, antisymmetric4(), i1, i2, i3, i4).hold();

	return indexed(pos_sig ? epsilon_pos : epsilon_neg, antisymmetric4(), i1, i2, i3, i4);
}

} // namespace GiNaC
