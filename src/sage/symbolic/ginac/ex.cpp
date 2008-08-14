/** @file ex.cpp
 *
 *  Implementation of GiNaC's light-weight expression handles. */

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

#include "ex.h"
#include "add.h"
#include "mul.h"
#include "ncmul.h"
#include "numeric.h"
#include "matrix.h"
#include "power.h"
#include "lst.h"
#include "relational.h"
#include "input_lexer.h"
#include "utils.h"

namespace GiNaC {

//////////
// other constructors
//////////

// none (all inlined)

//////////
// non-virtual functions in this class
//////////

// public
	
/** Print expression to stream. The formatting of the output is determined
 *  by the kind of print_context object that is passed. Possible formattings
 *  include ginsh-parsable output (the default), tree-like output for
 *  debugging, and C++ source.
 *  @see print_context */
void ex::print(const print_context & c, unsigned level) const
{
	bp->print(c, level);
}

/** Little wrapper arount print to be called within a debugger. */
void ex::dbgprint() const
{
	bp->dbgprint();
}

/** Little wrapper arount printtree to be called within a debugger. */
void ex::dbgprinttree() const
{
	bp->dbgprinttree();
}

ex ex::expand(unsigned options) const
{
	if (options == 0 && (bp->flags & status_flags::expanded)) // The "expanded" flag only covers the standard options; someone might want to re-expand with different options
		return *this;
	else
		return bp->expand(options);
}

/** Compute partial derivative of an expression.
 *
 *  @param s  symbol by which the expression is derived
 *  @param nth  order of derivative (default 1)
 *  @return partial derivative as a new expression */
ex ex::diff(const symbol & s, unsigned nth) const
{
	if (!nth)
		return *this;
	else
		return bp->diff(s, nth);
}

/** Check whether expression matches a specified pattern. */
bool ex::match(const ex & pattern) const
{
	lst repl_lst;
	return bp->match(pattern, repl_lst);
}

/** Find all occurrences of a pattern. The found matches are appended to
 *  the "found" list. If the expression itself matches the pattern, the
 *  children are not further examined. This function returns true when any
 *  matches were found. */
bool ex::find(const ex & pattern, lst & found) const
{
	if (match(pattern)) {
		found.append(*this);
		found.sort();
		found.unique();
		return true;
	}
	bool any_found = false;
	for (size_t i=0; i<nops(); i++)
		if (op(i).find(pattern, found))
			any_found = true;
	return any_found;
}

/** Substitute objects in an expression (syntactic substitution) and return
 *  the result as a new expression. */
ex ex::subs(const lst & ls, const lst & lr, unsigned options) const
{
	GINAC_ASSERT(ls.nops() == lr.nops());

	// Convert the lists to a map
	exmap m;
	for (lst::const_iterator its = ls.begin(), itr = lr.begin(); its != ls.end(); ++its, ++itr) {
		m.insert(std::make_pair(*its, *itr));

		// Search for products and powers in the expressions to be substituted
		// (for an optimization in expairseq::subs())
		if (is_exactly_a<mul>(*its) || is_exactly_a<power>(*its))
			options |= subs_options::pattern_is_product;
	}
	if (!(options & subs_options::pattern_is_product))
		options |= subs_options::pattern_is_not_product;

	return bp->subs(m, options);
}

/** Substitute objects in an expression (syntactic substitution) and return
 *  the result as a new expression.  There are two valid types of
 *  replacement arguments: 1) a relational like object==ex and 2) a list of
 *  relationals lst(object1==ex1,object2==ex2,...). */
ex ex::subs(const ex & e, unsigned options) const
{
	if (e.info(info_flags::relation_equal)) {

		// Argument is a relation: convert it to a map
		exmap m;
		const ex & s = e.op(0);
		m.insert(std::make_pair(s, e.op(1)));

		if (is_exactly_a<mul>(s) || is_exactly_a<power>(s))
			options |= subs_options::pattern_is_product;
		else
			options |= subs_options::pattern_is_not_product;

		return bp->subs(m, options);

	} else if (e.info(info_flags::list)) {

		// Argument is a list: convert it to a map
		exmap m;
		GINAC_ASSERT(is_a<lst>(e));
		for (lst::const_iterator it = ex_to<lst>(e).begin(); it != ex_to<lst>(e).end(); ++it) {
			ex r = *it;
			if (!r.info(info_flags::relation_equal))
				throw(std::invalid_argument("basic::subs(ex): argument must be a list of equations"));
			const ex & s = r.op(0);
			m.insert(std::make_pair(s, r.op(1)));

			// Search for products and powers in the expressions to be substituted
			// (for an optimization in expairseq::subs())
			if (is_exactly_a<mul>(s) || is_exactly_a<power>(s))
				options |= subs_options::pattern_is_product;
		}
		if (!(options & subs_options::pattern_is_product))
			options |= subs_options::pattern_is_not_product;

		return bp->subs(m, options);

	} else
		throw(std::invalid_argument("ex::subs(ex): argument must be a relation_equal or a list"));
}

/** Traverse expression tree with given visitor, preorder traversal. */
void ex::traverse_preorder(visitor & v) const
{
	accept(v);

	size_t n = nops();
	for (size_t i = 0; i < n; ++i)
		op(i).traverse_preorder(v);
}

/** Traverse expression tree with given visitor, postorder traversal. */
void ex::traverse_postorder(visitor & v) const
{
	size_t n = nops();
	for (size_t i = 0; i < n; ++i)
		op(i).traverse_postorder(v);

	accept(v);
}

/** Return modifyable operand/member at position i. */
ex & ex::let_op(size_t i)
{
	makewriteable();
	return bp->let_op(i);
}

ex & ex::operator[](const ex & index)
{
	makewriteable();
	return (*bp)[index];
}

ex & ex::operator[](size_t i)
{
	makewriteable();
	return (*bp)[i];
}

/** Left hand side of relational expression. */
ex ex::lhs() const
{
	if (!is_a<relational>(*this))
		throw std::runtime_error("ex::lhs(): not a relation");
	return bp->op(0);
}

/** Right hand side of relational expression. */
ex ex::rhs() const
{
	if (!is_a<relational>(*this))
		throw std::runtime_error("ex::rhs(): not a relation");
	return bp->op(1);
}

/** Check whether expression is a polynomial. */
bool ex::is_polynomial(const ex & vars) const
{
	if (is_a<lst>(vars)) {
		const lst & varlst = ex_to<lst>(vars);
		for (lst::const_iterator i=varlst.begin(); i!=varlst.end(); ++i)
			if (!bp->is_polynomial(*i))
				return false;
		return true;
	}
	else
		return bp->is_polynomial(vars);
}

/** Check whether expression is zero or zero matrix. */
bool ex::is_zero_matrix() const
{
	if (is_zero())
		return  true;
	else {
		ex e = evalm();
		return is_a<matrix>(e) && ex_to<matrix>(e).is_zero_matrix();
	}
}

// private

/** Make this ex writable (if more than one ex handle the same basic) by 
 *  unlinking the object and creating an unshared copy of it. */
void ex::makewriteable()
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
	bp.makewritable();
	GINAC_ASSERT(bp->get_refcount() == 1);
}

/** Share equal objects between expressions.
 *  @see ex::compare(const ex &) */
void ex::share(const ex & other) const
{
	if ((bp->flags | other.bp->flags) & status_flags::not_shareable)
		return;

	if (bp->get_refcount() <= other.bp->get_refcount())
		bp = other.bp;
	else
		other.bp = bp;
}

/** Helper function for the ex-from-basic constructor. This is where GiNaC's
 *  automatic evaluator and memory management are implemented.
 *  @see ex::ex(const basic &) */
ptr<basic> ex::construct_from_basic(const basic & other)
{
	if (!(other.flags & status_flags::evaluated)) {

		// The object is not yet evaluated, so call eval() to evaluate
		// the top level. This will return either
		//  a) the original object with status_flags::evaluated set (when the
		//     eval() implementation calls hold())
		// or
		//  b) a different expression.
		//
		// eval() returns an ex, not a basic&, so this will go through
		// construct_from_basic() a second time. In case a) we end up in
		// the "else" branch below. In case b) we end up here again and
		// apply eval() once more. The recursion stops when eval() calls
		// hold() or returns an object that already has its "evaluated"
		// flag set, such as a symbol or a numeric.
		const ex & tmpex = other.eval(1);

		// Eventually, the eval() recursion goes through the "else" branch
		// below, which assures that the object pointed to by tmpex.bp is
		// allocated on the heap (either it was already on the heap or it
		// is a heap-allocated duplicate of another object).
		GINAC_ASSERT(tmpex.bp->flags & status_flags::dynallocated); 

		// If the original object is not referenced but heap-allocated,
		// it means that eval() hit case b) above. The original object is
		// no longer needed (it evaluated into something different), so we
		// delete it (because nobody else will).
		if ((other.get_refcount() == 0) && (other.flags & status_flags::dynallocated))
			delete &other; // yes, you can apply delete to a const pointer

		// We can't return a basic& here because the tmpex is destroyed as
		// soon as we leave the function, which would deallocate the
		// evaluated object.
		return tmpex.bp;

	} else {

		// The easy case: making an "ex" out of an evaluated object.
		if (other.flags & status_flags::dynallocated) {

			// The object is already heap-allocated, so we can just make
			// another reference to it.
			return ptr<basic>(const_cast<basic &>(other));

		} else {

			// The object is not heap-allocated, so we create a duplicate
			// on the heap.
			basic *bp = other.duplicate();
			bp->setflag(status_flags::dynallocated);
			GINAC_ASSERT(bp->get_refcount() == 0);
			return bp;
		}
	}
}

basic & ex::construct_from_int(int i)
{
	switch (i) {  // prefer flyweights over new objects
	case -12:
		return *const_cast<numeric *>(_num_12_p);
	case -11:
		return *const_cast<numeric *>(_num_11_p);
	case -10:
		return *const_cast<numeric *>(_num_10_p);
	case -9:
		return *const_cast<numeric *>(_num_9_p);
	case -8:
		return *const_cast<numeric *>(_num_8_p);
	case -7:
		return *const_cast<numeric *>(_num_7_p);
	case -6:
		return *const_cast<numeric *>(_num_6_p);
	case -5:
		return *const_cast<numeric *>(_num_5_p);
	case -4:
		return *const_cast<numeric *>(_num_4_p);
	case -3:
		return *const_cast<numeric *>(_num_3_p);
	case -2:
		return *const_cast<numeric *>(_num_2_p);
	case -1:
		return *const_cast<numeric *>(_num_1_p);
	case 0:
		return *const_cast<numeric *>(_num0_p);
	case 1:
		return *const_cast<numeric *>(_num1_p);
	case 2:
		return *const_cast<numeric *>(_num2_p);
	case 3:
		return *const_cast<numeric *>(_num3_p);
	case 4:
		return *const_cast<numeric *>(_num4_p);
	case 5:
		return *const_cast<numeric *>(_num5_p);
	case 6:
		return *const_cast<numeric *>(_num6_p);
	case 7:
		return *const_cast<numeric *>(_num7_p);
	case 8:
		return *const_cast<numeric *>(_num8_p);
	case 9:
		return *const_cast<numeric *>(_num9_p);
	case 10:
		return *const_cast<numeric *>(_num10_p);
	case 11:
		return *const_cast<numeric *>(_num11_p);
	case 12:
		return *const_cast<numeric *>(_num12_p);
	default:
		basic *bp = new numeric(i);
		bp->setflag(status_flags::dynallocated);
		GINAC_ASSERT(bp->get_refcount() == 0);
		return *bp;
	}
}
	
basic & ex::construct_from_uint(unsigned int i)
{
	switch (i) {  // prefer flyweights over new objects
	case 0:
		return *const_cast<numeric *>(_num0_p);
	case 1:
		return *const_cast<numeric *>(_num1_p);
	case 2:
		return *const_cast<numeric *>(_num2_p);
	case 3:
		return *const_cast<numeric *>(_num3_p);
	case 4:
		return *const_cast<numeric *>(_num4_p);
	case 5:
		return *const_cast<numeric *>(_num5_p);
	case 6:
		return *const_cast<numeric *>(_num6_p);
	case 7:
		return *const_cast<numeric *>(_num7_p);
	case 8:
		return *const_cast<numeric *>(_num8_p);
	case 9:
		return *const_cast<numeric *>(_num9_p);
	case 10:
		return *const_cast<numeric *>(_num10_p);
	case 11:
		return *const_cast<numeric *>(_num11_p);
	case 12:
		return *const_cast<numeric *>(_num12_p);
	default:
		basic *bp = new numeric(i);
		bp->setflag(status_flags::dynallocated);
		GINAC_ASSERT(bp->get_refcount() == 0);
		return *bp;
	}
}
	
basic & ex::construct_from_long(long i)
{
	switch (i) {  // prefer flyweights over new objects
	case -12:
		return *const_cast<numeric *>(_num_12_p);
	case -11:
		return *const_cast<numeric *>(_num_11_p);
	case -10:
		return *const_cast<numeric *>(_num_10_p);
	case -9:
		return *const_cast<numeric *>(_num_9_p);
	case -8:
		return *const_cast<numeric *>(_num_8_p);
	case -7:
		return *const_cast<numeric *>(_num_7_p);
	case -6:
		return *const_cast<numeric *>(_num_6_p);
	case -5:
		return *const_cast<numeric *>(_num_5_p);
	case -4:
		return *const_cast<numeric *>(_num_4_p);
	case -3:
		return *const_cast<numeric *>(_num_3_p);
	case -2:
		return *const_cast<numeric *>(_num_2_p);
	case -1:
		return *const_cast<numeric *>(_num_1_p);
	case 0:
		return *const_cast<numeric *>(_num0_p);
	case 1:
		return *const_cast<numeric *>(_num1_p);
	case 2:
		return *const_cast<numeric *>(_num2_p);
	case 3:
		return *const_cast<numeric *>(_num3_p);
	case 4:
		return *const_cast<numeric *>(_num4_p);
	case 5:
		return *const_cast<numeric *>(_num5_p);
	case 6:
		return *const_cast<numeric *>(_num6_p);
	case 7:
		return *const_cast<numeric *>(_num7_p);
	case 8:
		return *const_cast<numeric *>(_num8_p);
	case 9:
		return *const_cast<numeric *>(_num9_p);
	case 10:
		return *const_cast<numeric *>(_num10_p);
	case 11:
		return *const_cast<numeric *>(_num11_p);
	case 12:
		return *const_cast<numeric *>(_num12_p);
	default:
		basic *bp = new numeric(i);
		bp->setflag(status_flags::dynallocated);
		GINAC_ASSERT(bp->get_refcount() == 0);
		return *bp;
	}
}
	
basic & ex::construct_from_ulong(unsigned long i)
{
	switch (i) {  // prefer flyweights over new objects
	case 0:
		return *const_cast<numeric *>(_num0_p);
	case 1:
		return *const_cast<numeric *>(_num1_p);
	case 2:
		return *const_cast<numeric *>(_num2_p);
	case 3:
		return *const_cast<numeric *>(_num3_p);
	case 4:
		return *const_cast<numeric *>(_num4_p);
	case 5:
		return *const_cast<numeric *>(_num5_p);
	case 6:
		return *const_cast<numeric *>(_num6_p);
	case 7:
		return *const_cast<numeric *>(_num7_p);
	case 8:
		return *const_cast<numeric *>(_num8_p);
	case 9:
		return *const_cast<numeric *>(_num9_p);
	case 10:
		return *const_cast<numeric *>(_num10_p);
	case 11:
		return *const_cast<numeric *>(_num11_p);
	case 12:
		return *const_cast<numeric *>(_num12_p);
	default:
		basic *bp = new numeric(i);
		bp->setflag(status_flags::dynallocated);
		GINAC_ASSERT(bp->get_refcount() == 0);
		return *bp;
	}
}
	
basic & ex::construct_from_double(double d)
{
	basic *bp = new numeric(d);
	bp->setflag(status_flags::dynallocated);
	GINAC_ASSERT(bp->get_refcount() == 0);
	return *bp;
}

ptr<basic> ex::construct_from_string_and_lst(const std::string &s, const ex &l)
{
	set_lexer_string(s);
	set_lexer_symbols(l);
	ginac_yyrestart(NULL);
	if (ginac_yyparse())
		throw (std::runtime_error(get_parser_error()));
	else
		return parsed_ex.bp;
}
	
//////////
// static member variables
//////////

// none

//////////
// functions which are not member functions
//////////

// none

//////////
// global functions
//////////

// none


} // namespace GiNaC
