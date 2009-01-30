/** @file power.cpp
 *
 *  Implementation of GiNaC's symbolic exponentiation (basis^exponent). */

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


// ADDED FOR SAGE; this gets us nearly 20% speedup already for (x+y+z)^10.
#include "Python.h"
extern "C" PyObject* py_binomial_int(int n, unsigned int k);


#include <vector>
#include <iostream>
#include <stdexcept>
#include <limits>

#include "power.h"
#include "expairseq.h"
#include "add.h"
#include "mul.h"
#include "ncmul.h"
#include "numeric.h"
#include "constant.h"
#include "operators.h"
#include "inifcns.h" // for log() in power::derivative()
#include "matrix.h"
#include "indexed.h"
#include "symbol.h"
#include "lst.h"
#include "archive.h"
#include "utils.h"
#include "relational.h"
#include "compiler.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(power, basic,
  print_func<print_dflt>(&power::do_print_dflt).
  print_func<print_latex>(&power::do_print_latex).
  print_func<print_csrc>(&power::do_print_csrc).
  print_func<print_python>(&power::do_print_python).
  print_func<print_python_repr>(&power::do_print_python_repr))

typedef std::vector<int> intvector;

//////////
// default constructor
//////////

power::power() : inherited(&power::tinfo_static) { }

//////////
// other constructors
//////////

// all inlined

//////////
// archiving
//////////

power::power(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	n.find_ex("basis", basis, sym_lst);
	n.find_ex("exponent", exponent, sym_lst);
}

void power::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_ex("basis", basis);
	n.add_ex("exponent", exponent);
}

DEFAULT_UNARCHIVE(power)

//////////
// functions overriding virtual functions from base classes
//////////

// public

unsigned power::calchash() const
{
  // TODO: SAGE -- we have to disable hashing in order to
  // get internal sorting of poly vars (via PARI) 
  // to work correctly on the linear part.  This makes
  // the hashes fail, so the actual compare funtion
  // is called.
  hashvalue = 0;
  setflag(status_flags::hash_calculated);
  return hashvalue;
}

void power::print_power(const print_context & c, const char *powersymbol, const char *openbrace, const char *closebrace, unsigned level) const
{
	// Ordinary output of powers using '^' or '**'
	if (precedence() <= level)
		c.s << openbrace << '(';
	basis.print(c, precedence());
	c.s << powersymbol;
	c.s << openbrace;
	exponent.print(c, precedence());
	c.s << closebrace;
	if (precedence() <= level)
		c.s << ')' << closebrace;
}

void power::do_print_dflt(const print_dflt & c, unsigned level) const
{
	if (exponent.is_equal(_ex1_2)) {

		// Square roots are printed in a special way
		c.s << "sqrt(";
		basis.print(c);
		c.s << ')';

	} else
		print_power(c, "^", "", "", level);
}

void power::do_print_latex(const print_latex & c, unsigned level) const
{
	if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_negative()) {

		// Powers with negative numeric exponents are printed as fractions
		c.s << "\\frac{1}{";
		power(basis, -exponent).eval().print(c);
		c.s << '}';

	} else if (exponent.is_equal(_ex1_2)) {

		// Square roots are printed in a special way
		c.s << "\\sqrt{";
		basis.print(c);
		c.s << '}';

	} else
		print_power(c, "^", "{", "}", level);
}

static void print_sym_pow(const print_context & c, const symbol &x, int exp)
{
	// Optimal output of integer powers of symbols to aid compiler CSE.
	// C.f. ISO/IEC 14882:1998, section 1.9 [intro execution], paragraph 15
	// to learn why such a parenthesation is really necessary.
	if (exp == 1) {
		x.print(c);
	} else if (exp == 2) {
		x.print(c);
		c.s << "*";
		x.print(c);
	} else if (exp & 1) {
		x.print(c);
		c.s << "*";
		print_sym_pow(c, x, exp-1);
	} else {
		c.s << "(";
		print_sym_pow(c, x, exp >> 1);
		c.s << ")*(";
		print_sym_pow(c, x, exp >> 1);
		c.s << ")";
	}
}

void power::do_print_csrc(const print_csrc & c, unsigned level) const
{
	if (is_a<print_csrc_cl_N>(c)) {
		if (exponent.is_equal(_ex_1)) {
			c.s << "recip(";
			basis.print(c);
			c.s << ')';
			return;
		}
		c.s << "expt(";
		basis.print(c);
		c.s << ", ";
		exponent.print(c);
		c.s << ')';
		return;
	}

	// Integer powers of symbols are printed in a special, optimized way
	if (exponent.info(info_flags::integer)
	 && (is_a<symbol>(basis) || is_a<constant>(basis))) {
		int exp = ex_to<numeric>(exponent).to_int();
		if (exp > 0)
			c.s << '(';
		else {
			exp = -exp;
			c.s << "1.0/(";
		}
		print_sym_pow(c, ex_to<symbol>(basis), exp);
		c.s << ')';

	// <expr>^-1 is printed as "1.0/<expr>" or with the recip() function of CLN
	} else if (exponent.is_equal(_ex_1)) {
		c.s << "1.0/(";
		basis.print(c);
		c.s << ')';

	// Otherwise, use the pow() function
	} else {
		c.s << "pow(";
		basis.print(c);
		c.s << ',';
		exponent.print(c);
		c.s << ')';
	}
}

void power::do_print_python(const print_python & c, unsigned level) const
{
	print_power(c, "**", "", "", level);
}

void power::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(';
	basis.print(c);
	c.s << ',';
	exponent.print(c);
	c.s << ')';
}

bool power::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::polynomial:
		case info_flags::integer_polynomial:
		case info_flags::cinteger_polynomial:
		case info_flags::rational_polynomial:
		case info_flags::crational_polynomial:
			return exponent.info(info_flags::nonnegint) &&
			       basis.info(inf);
		case info_flags::rational_function:
			return exponent.info(info_flags::integer) &&
			       basis.info(inf);
		case info_flags::algebraic:
			return !exponent.info(info_flags::integer) ||
			       basis.info(inf);
		case info_flags::expanded:
			return (flags & status_flags::expanded);
		case info_flags::has_indices: {
			if (flags & status_flags::has_indices)
				return true;
			else if (flags & status_flags::has_no_indices)
				return false;
			else if (basis.info(info_flags::has_indices)) {
				setflag(status_flags::has_indices);
				clearflag(status_flags::has_no_indices);
				return true;
			} else {
				clearflag(status_flags::has_indices);
				setflag(status_flags::has_no_indices);
				return false;
			}
		}
	}
	return inherited::info(inf);
}

size_t power::nops() const
{
	return 2;
}

ex power::op(size_t i) const
{
	GINAC_ASSERT(i<2);

	return i==0 ? basis : exponent;
}

ex power::map(map_function & f) const
{
	const ex &mapped_basis = f(basis);
	const ex &mapped_exponent = f(exponent);

	if (!are_ex_trivially_equal(basis, mapped_basis)
	 || !are_ex_trivially_equal(exponent, mapped_exponent))
		return (new power(mapped_basis, mapped_exponent))->setflag(status_flags::dynallocated);
	else
		return *this;
}

bool power::is_polynomial(const ex & var) const
{
	if (exponent.has(var))
		return false;
	if (!exponent.info(info_flags::nonnegint))
		return false;
	return basis.is_polynomial(var);
}

int power::degree(const ex & s) const
{
	if (is_equal(ex_to<basic>(s)))
		return 1;
	else if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_integer()) {
		if (basis.is_equal(s))
			return ex_to<numeric>(exponent).to_int();
		else
			return basis.degree(s) * ex_to<numeric>(exponent).to_int();
	} else if (basis.has(s))
		throw(std::runtime_error("power::degree(): undefined degree because of non-integer exponent"));
	else
		return 0;
}

int power::ldegree(const ex & s) const 
{
	if (is_equal(ex_to<basic>(s)))
		return 1;
	else if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_integer()) {
		if (basis.is_equal(s))
			return ex_to<numeric>(exponent).to_int();
		else
			return basis.ldegree(s) * ex_to<numeric>(exponent).to_int();
	} else if (basis.has(s))
		throw(std::runtime_error("power::ldegree(): undefined degree because of non-integer exponent"));
	else
		return 0;
}

ex power::coeff(const ex & s, int n) const
{
	if (is_equal(ex_to<basic>(s)))
		return n==1 ? _ex1 : _ex0;
	else if (!basis.is_equal(s)) {
		// basis not equal to s
		if (n == 0)
			return *this;
		else
			return _ex0;
	} else {
		// basis equal to s
		if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_integer()) {
			// integer exponent
			int int_exp = ex_to<numeric>(exponent).to_int();
			if (n == int_exp)
				return _ex1;
			else
				return _ex0;
		} else {
			// non-integer exponents are treated as zero
			if (n == 0)
				return *this;
			else
				return _ex0;
		}
	}
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x, x1, x2,... stand for a symbolic variables of type ex and c, c1, c2...
 *  stand for such expressions that contain a plain number.
 *  - ^(x,0) -> 1  (also handles ^(0,0))
 *  - ^(x,1) -> x
 *  - ^(0,c) -> 0 or exception  (depending on the real part of c)
 *  - ^(1,x) -> 1
 *  - ^(c1,c2) -> *(c1^n,c1^(c2-n))  (so that 0<(c2-n)<1, try to evaluate roots, possibly in numerator and denominator of c1)
 *  - ^(^(x,c1),c2) -> ^(x,c1*c2)  if x is positive and c1 is real.
 *  - ^(^(x,c1),c2) -> ^(x,c1*c2)  (c2 integer or -1 < c1 <= 1, case c1=1 should not happen, see below!)
 *  - ^(*(x,y,z),c) -> *(x^c,y^c,z^c)  (if c integer)
 *  - ^(*(x,c1),c2) -> ^(x,c2)*c1^c2  (c1>0)
 *  - ^(*(x,c1),c2) -> ^(-x,c2)*c1^c2  (c1<0)
 *
 *  @param level cut-off in recursive evaluation */
ex power::eval(int level) const
{
	if ((level==1) && (flags & status_flags::evaluated))
		return *this;
	else if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	const ex & ebasis    = level==1 ? basis    : basis.eval(level-1);
	const ex & eexponent = level==1 ? exponent : exponent.eval(level-1);
	
	bool basis_is_numerical = false;
	bool exponent_is_numerical = false;
	const numeric *num_basis;
	const numeric *num_exponent;
	
	if (is_exactly_a<numeric>(ebasis)) {
		basis_is_numerical = true;
		num_basis = &ex_to<numeric>(ebasis);
	}
	if (is_exactly_a<numeric>(eexponent)) {
		exponent_is_numerical = true;
		num_exponent = &ex_to<numeric>(eexponent);
	}
	
	// ^(x,0) -> 1  (0^0 also handled here)
	if (eexponent.is_zero() && 
		!(is_exactly_a<numeric>(ebasis) && 
			ex_to<numeric>(ebasis).is_parent_pos_char())) {
		if (ebasis.is_zero())
			throw (std::domain_error("power::eval(): pow(0,0) is undefined"));
		else
			return _ex1;
	}
	
	// ^(x,1) -> x
	if (eexponent.is_equal(_ex1))
		return ebasis;

	// ^(0,c1) -> 0 or exception  (depending on real value of c1)
	if (ebasis.is_zero() && exponent_is_numerical) {
		if ((num_exponent->real()).is_zero())
			throw (std::domain_error("power::eval(): pow(0,I) is undefined"));
		else if ((num_exponent->real()).is_negative())
			throw (pole_error("power::eval(): division by zero",1));
		else
			return _ex0;
	}

	// ^(1,x) -> 1
	if (ebasis.is_equal(_ex1))
		return _ex1;

	// power of a function calculated by separate rules defined for this function
	if (is_exactly_a<function>(ebasis))
		return ex_to<function>(ebasis).power(eexponent);

	// Turn (x^c)^d into x^(c*d) in the case that x is positive and c is real.
	if (is_exactly_a<power>(ebasis) && ebasis.op(0).info(info_flags::positive) && ebasis.op(1).info(info_flags::real))
		return power(ebasis.op(0), ebasis.op(1) * eexponent);

	if (exponent_is_numerical) {

		// ^(c1,c2) -> c1^c2  (c1, c2 numeric(),
		// except if c1,c2 are rational, but c1^c2 is not)
		if (basis_is_numerical) {
			const bool basis_is_crational = num_basis->is_crational();
			const bool exponent_is_crational = num_exponent->is_crational();
			if (!basis_is_crational || !exponent_is_crational) {
				// return a plain float
				return (new numeric(num_basis->power(*num_exponent)))->setflag(status_flags::dynallocated |
				                                                               status_flags::evaluated |
				                                                               status_flags::expanded);
			}

			const numeric res = num_basis->power(*num_exponent);
			if (res.is_crational()) {
				return res;
			}
			GINAC_ASSERT(!num_exponent->is_integer());  // has been handled by now

			// ^(c1,n/m) -> *(c1^q,c1^(n/m-q)), 0<(n/m-q)<1, q integer
			if (basis_is_crational && exponent_is_crational
			    && num_exponent->is_real()
			    && !num_exponent->is_integer()) {
				const numeric n = num_exponent->numer();
				const numeric m = num_exponent->denom();
				numeric r;
				numeric q = iquo(n, m, r);
				if (r.is_negative()) {
					r += m;
					--q;
				}
				if (q.is_zero()) {  // the exponent was in the allowed range 0<(n/m)<1
					if (num_basis->is_rational() && !num_basis->is_integer()) {
						// try it for numerator and denominator separately, in order to
						// partially simplify things like (5/8)^(1/3) -> 1/2*5^(1/3)
						const numeric bnum = num_basis->numer();
						const numeric bden = num_basis->denom();
						const numeric res_bnum = bnum.power(*num_exponent);
						const numeric res_bden = bden.power(*num_exponent);
						if (res_bnum.is_integer())
							return (new mul(power(bden,-*num_exponent),res_bnum))->setflag(status_flags::dynallocated | status_flags::evaluated);
						if (res_bden.is_integer())
							return (new mul(power(bnum,*num_exponent),res_bden.inverse()))->setflag(status_flags::dynallocated | status_flags::evaluated);
					}
					return this->hold();
				} else {
					// assemble resulting product, but allowing for a re-evaluation,
					// because otherwise we'll end up with something like
					//    (7/8)^(4/3)  ->  7/8*(1/2*7^(1/3))
					// instead of 7/16*7^(1/3).
					ex prod = power(*num_basis,r.div(m));
					return prod*power(*num_basis,q);
				}
			}
		}
	
		// ^(^(x,c1),c2) -> ^(x,c1*c2)
		// (c1, c2 numeric(), c2 integer or -1 < c1 <= 1,
		// case c1==1 should not happen, see below!)
		if (is_exactly_a<power>(ebasis)) {
			const power & sub_power = ex_to<power>(ebasis);
			const ex & sub_basis = sub_power.basis;
			const ex & sub_exponent = sub_power.exponent;
			if (is_exactly_a<numeric>(sub_exponent)) {
				const numeric & num_sub_exponent = ex_to<numeric>(sub_exponent);
				GINAC_ASSERT(num_sub_exponent!=numeric(1));
				if (num_exponent->is_integer() || (abs(num_sub_exponent) - (*_num1_p)).is_negative()) {
					return power(sub_basis,num_sub_exponent.mul(*num_exponent));
				}
			}
		}
	
		// ^(*(x,y,z),c1) -> *(x^c1,y^c1,z^c1) (c1 integer)
		if (num_exponent->is_integer() && is_exactly_a<mul>(ebasis)) {
			return expand_mul(ex_to<mul>(ebasis), *num_exponent, 0);
		}

		// (2*x + 6*y)^(-4) -> 1/16*(x + 3*y)^(-4)
		if (num_exponent->is_integer() && is_exactly_a<add>(ebasis)) {
			numeric icont = ebasis.integer_content();
			const numeric lead_coeff = 
				ex_to<numeric>(ex_to<add>(ebasis).seq.begin()->coeff).div(icont);

			const bool canonicalizable = lead_coeff.is_integer();
			const bool unit_normal = lead_coeff.is_pos_integer();
			if (canonicalizable && (! unit_normal))
				icont = icont.mul(*_num_1_p);
			
			if (canonicalizable && (icont != *_num1_p)) {
				const add& addref = ex_to<add>(ebasis);
				add* addp = new add(addref);
				addp->setflag(status_flags::dynallocated);
				addp->clearflag(status_flags::hash_calculated);
				addp->overall_coeff = ex_to<numeric>(addp->overall_coeff).div_dyn(icont);
				for (epvector::iterator i = addp->seq.begin(); i != addp->seq.end(); ++i)
					i->coeff = ex_to<numeric>(i->coeff).div_dyn(icont);

				const numeric c = icont.power(*num_exponent);
				if (likely(c != *_num1_p))
					return (new mul(power(*addp, *num_exponent), c))->setflag(status_flags::dynallocated);
				else
					return power(*addp, *num_exponent);
			}
		}

		// ^(*(...,x;c1),c2) -> *(^(*(...,x;1),c2),c1^c2)  (c1, c2 numeric(), c1>0)
		// ^(*(...,x;c1),c2) -> *(^(*(...,x;-1),c2),(-c1)^c2)  (c1, c2 numeric(), c1<0)
		if (is_exactly_a<mul>(ebasis)) {
			GINAC_ASSERT(!num_exponent->is_integer()); // should have been handled above
			const mul & mulref = ex_to<mul>(ebasis);
			if (!mulref.overall_coeff.is_equal(_ex1)) {
				const numeric & num_coeff = ex_to<numeric>(mulref.overall_coeff);
				if (num_coeff.is_real()) {
					if (num_coeff.is_positive()) {
						mul *mulp = new mul(mulref);
						mulp->overall_coeff = _ex1;
						mulp->clearflag(status_flags::evaluated);
						mulp->clearflag(status_flags::hash_calculated);
						return (new mul(power(*mulp,exponent),
						                power(num_coeff,*num_exponent)))->setflag(status_flags::dynallocated);
					} else {
						GINAC_ASSERT(num_coeff.compare(*_num0_p)<0);
						if (!num_coeff.is_equal(*_num_1_p)) {
							mul *mulp = new mul(mulref);
							mulp->overall_coeff = _ex_1;
							mulp->clearflag(status_flags::evaluated);
							mulp->clearflag(status_flags::hash_calculated);
							return (new mul(power(*mulp,exponent),
							                power(abs(num_coeff),*num_exponent)))->setflag(status_flags::dynallocated);
						}
					}
				}
			}
		}

		// ^(nc,c1) -> ncmul(nc,nc,...) (c1 positive integer, unless nc is a matrix)
		if (num_exponent->is_pos_integer() &&
		    ebasis.return_type() != return_types::commutative &&
		    !is_a<matrix>(ebasis)) {
			return ncmul(exvector(num_exponent->to_int(), ebasis), true);
		}
	}
	
	if (are_ex_trivially_equal(ebasis,basis) &&
	    are_ex_trivially_equal(eexponent,exponent)) {
		return this->hold();
	}
	return (new power(ebasis, eexponent))->setflag(status_flags::dynallocated |
	                                               status_flags::evaluated);
}

ex power::evalf(int level) const
{
	ex ebasis;
	ex eexponent;
	
	if (level==1) {
		ebasis = basis;
		eexponent = exponent;
	} else if (level == -max_recursion_level) {
		throw(std::runtime_error("max recursion level reached"));
	} else {
		ebasis = basis.evalf(level-1);
		if (!is_exactly_a<numeric>(exponent))
			eexponent = exponent.evalf(level-1);
		else
			eexponent = exponent;
	}

	return power(ebasis,eexponent);
}

ex power::evalm() const
{
	const ex ebasis = basis.evalm();
	const ex eexponent = exponent.evalm();
	if (is_a<matrix>(ebasis)) {
		if (is_exactly_a<numeric>(eexponent)) {
			return (new matrix(ex_to<matrix>(ebasis).pow(eexponent)))->setflag(status_flags::dynallocated);
		}
	}
	return (new power(ebasis, eexponent))->setflag(status_flags::dynallocated);
}

bool power::has(const ex & other, unsigned options) const
{
	if (!(options & has_options::algebraic))
		return basic::has(other, options);
	if (!is_a<power>(other))
		return basic::has(other, options);
	if (!exponent.info(info_flags::integer)
			|| !other.op(1).info(info_flags::integer))
		return basic::has(other, options);
	if (exponent.info(info_flags::posint)
			&& other.op(1).info(info_flags::posint)
			&& ex_to<numeric>(exponent).to_int()
					> ex_to<numeric>(other.op(1)).to_int()
			&& basis.match(other.op(0)))
		return true;
	if (exponent.info(info_flags::negint)
			&& other.op(1).info(info_flags::negint)
			&& ex_to<numeric>(exponent).to_int()
					< ex_to<numeric>(other.op(1)).to_int()
			&& basis.match(other.op(0)))
		return true;
	return basic::has(other, options);
}

// from mul.cpp
extern bool tryfactsubs(const ex &, const ex &, int &, lst &);

ex power::subs(const exmap & m, unsigned options) const
{	
	const ex &subsed_basis = basis.subs(m, options);
	const ex &subsed_exponent = exponent.subs(m, options);

	if (!are_ex_trivially_equal(basis, subsed_basis)
	 || !are_ex_trivially_equal(exponent, subsed_exponent)) 
		return power(subsed_basis, subsed_exponent).subs_one_level(m, options);

	if (!(options & subs_options::algebraic))
		return subs_one_level(m, options);

	for (exmap::const_iterator it = m.begin(); it != m.end(); ++it) {
		int nummatches = std::numeric_limits<int>::max();
		lst repls;
		if (tryfactsubs(*this, it->first, nummatches, repls))
			return (ex_to<basic>((*this) * power(it->second.subs(ex(repls), subs_options::no_pattern) / it->first.subs(ex(repls), subs_options::no_pattern), nummatches))).subs_one_level(m, options);
	}

	return subs_one_level(m, options);
}

ex power::eval_ncmul(const exvector & v) const
{
	return inherited::eval_ncmul(v);
}

ex power::conjugate() const
{
	ex newbasis = basis.conjugate();
	ex newexponent = exponent.conjugate();
	if (are_ex_trivially_equal(basis, newbasis) && are_ex_trivially_equal(exponent, newexponent)) {
		return *this;
	}
	return (new power(newbasis, newexponent))->setflag(status_flags::dynallocated);
}

ex power::real_part() const
{
	if (exponent.info(info_flags::integer)) {
		ex basis_real = basis.real_part();
		if (basis_real == basis)
			return *this;
		realsymbol a("a"),b("b");
		ex result;
		if (exponent.info(info_flags::posint))
			result = power(a+I*b,exponent);
		else
			result = power(a/(a*a+b*b)-I*b/(a*a+b*b),-exponent);
		result = result.expand();
		result = result.real_part();
		result = result.subs(lst( a==basis_real, b==basis.imag_part() ));
		return result;
	}
	
	ex a = basis.real_part();
	ex b = basis.imag_part();
	ex c = exponent.real_part();
	ex d = exponent.imag_part();
	return power(abs(basis),c)*exp(-d*atan2(b,a))*cos(c*atan2(b,a)+d*log(abs(basis)));
}

ex power::imag_part() const
{
	if (exponent.info(info_flags::integer)) {
		ex basis_real = basis.real_part();
		if (basis_real == basis)
			return 0;
		realsymbol a("a"),b("b");
		ex result;
		if (exponent.info(info_flags::posint))
			result = power(a+I*b,exponent);
		else
			result = power(a/(a*a+b*b)-I*b/(a*a+b*b),-exponent);
		result = result.expand();
		result = result.imag_part();
		result = result.subs(lst( a==basis_real, b==basis.imag_part() ));
		return result;
	}
	
	ex a=basis.real_part();
	ex b=basis.imag_part();
	ex c=exponent.real_part();
	ex d=exponent.imag_part();
	return
		power(abs(basis),c)*exp(-d*atan2(b,a))*sin(c*atan2(b,a)+d*log(abs(basis)));
}

// protected

// protected

/** Implementation of ex::diff() for a power.
 *  @see ex::diff */
ex power::derivative(const symbol & s) const
{
	if (is_a<numeric>(exponent)) {
		// D(b^r) = r * b^(r-1) * D(b) (faster than the formula below)
		epvector newseq;
		newseq.reserve(2);
		newseq.push_back(expair(basis, exponent - _ex1));
		newseq.push_back(expair(basis.diff(s), _ex1));
		return mul(newseq, exponent);
	} else {
		// D(b^e) = b^e * (D(e)*ln(b) + e*D(b)/b)
		return mul(*this,
		           add(mul(exponent.diff(s), log(basis)),
		           mul(mul(exponent, basis.diff(s)), power(basis, _ex_1))));
	}
}

int power::compare(const basic& other) const
{
	static const tinfo_t mul_id = find_tinfo_key("mul");
	static const tinfo_t symbol_id = find_tinfo_key("symbol");
	static const tinfo_t function_id = find_tinfo_key("function");
	static const tinfo_t fderivative_id = find_tinfo_key("fderivative");
	const tinfo_t typeid_this = tinfo();
	const tinfo_t typeid_other = other.tinfo();
	if (typeid_this==typeid_other) {
		GINAC_ASSERT(typeid(*this)==typeid(other));
		return compare_same_type(other);
	} else if (typeid_other == mul_id) {
		return -static_cast<const mul&>(other).compare_pow(*this);
	} else if (typeid_other == symbol_id) {
		return compare_symbol(static_cast<const symbol&>(other));
	} else if (typeid_other == function_id ||
			typeid_other == fderivative_id) {
		return -1;
	} else {
		return (typeid_this<typeid_other ? -1 : 1);
	}
}

int power::compare_symbol(const symbol & other) const
{
	int cmpval;
	cmpval = _ex1.compare(exponent);
	if (cmpval != 0) {
		return cmpval;
	}
	return basis.compare(other);
}

int power::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<power>(other));
	const power &o = static_cast<const power &>(other);

	int cmpval = basis.compare(o.basis);
	if (cmpval)
		return cmpval;
	else
	  // SAGE -- I changed the sign below for consistency with standard
	  // mathematics and all other math software (except mathematica).
	  // This makes it so x^5 + x^2 prints correctly instead of 
	  // as x^2 + x^5.   -- William Stein
		return -exponent.compare(o.exponent);
}

unsigned power::return_type() const
{
	return basis.return_type();
}

tinfo_t power::return_type_tinfo() const
{
	return basis.return_type_tinfo();
}

ex power::expand(unsigned options) const
{
	if (is_a<symbol>(basis) && exponent.info(info_flags::integer)) {
		// A special case worth optimizing.
		setflag(status_flags::expanded);
		return *this;
	}

	const ex expanded_basis = basis.expand(options);
	const ex expanded_exponent = exponent.expand(options);
	
	// x^(a+b) -> x^a * x^b
	if (is_exactly_a<add>(expanded_exponent)) {
		const add &a = ex_to<add>(expanded_exponent);
		exvector distrseq;
		distrseq.reserve(a.seq.size() + 1);
		epvector::const_iterator last = a.seq.end();
		epvector::const_iterator cit = a.seq.begin();
		while (cit!=last) {
			distrseq.push_back(power(expanded_basis, a.recombine_pair_to_ex(*cit)));
			++cit;
		}
		
		// Make sure that e.g. (x+y)^(2+a) expands the (x+y)^2 factor
		if (ex_to<numeric>(a.overall_coeff).is_integer()) {
			const numeric &num_exponent = ex_to<numeric>(a.overall_coeff);
			int int_exponent = num_exponent.to_int();
			if (int_exponent > 0 && is_exactly_a<add>(expanded_basis))
				distrseq.push_back(expand_add(ex_to<add>(expanded_basis), int_exponent, options));
			else
				distrseq.push_back(power(expanded_basis, a.overall_coeff));
		} else
			distrseq.push_back(power(expanded_basis, a.overall_coeff));
		
		// Make sure that e.g. (x+y)^(1+a) -> x*(x+y)^a + y*(x+y)^a
		ex r = (new mul(distrseq))->setflag(status_flags::dynallocated);
		return r.expand(options);
	}
	
	if (!is_exactly_a<numeric>(expanded_exponent) ||
		!ex_to<numeric>(expanded_exponent).is_integer()) {
		if (are_ex_trivially_equal(basis,expanded_basis) && are_ex_trivially_equal(exponent,expanded_exponent)) {
			return this->hold();
		} else {
			return (new power(expanded_basis,expanded_exponent))->setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
		}
	}
	
	// integer numeric exponent
	const numeric & num_exponent = ex_to<numeric>(expanded_exponent);
	int int_exponent = num_exponent.to_int();
	
	// (x+y)^n, n>0
	if (int_exponent > 0 && is_exactly_a<add>(expanded_basis))
		return expand_add(ex_to<add>(expanded_basis), int_exponent, options);
	
	// (x*y)^n -> x^n * y^n
	if (is_exactly_a<mul>(expanded_basis))
		return expand_mul(ex_to<mul>(expanded_basis), num_exponent, options, true);
	
	// cannot expand further
	if (are_ex_trivially_equal(basis,expanded_basis) && are_ex_trivially_equal(exponent,expanded_exponent))
		return this->hold();
	else
		return (new power(expanded_basis,expanded_exponent))->setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

/** expand a^n where a is an add and n is a positive integer.
 *  @see power::expand */
ex power::expand_add(const add & a, int n, unsigned options) const
{
	if (n==2)
		return expand_add_2(a, options);

	const size_t m = a.nops();
	exvector result;
	// The number of terms will be the number of combinatorial compositions,
	// i.e. the number of unordered arrangements of m nonnegative integers
	// which sum up to n.  It is frequently written as C_n(m) and directly
	// related with binomial coefficients:
	result.reserve(numeric(py_binomial_int(n+m-1, m-1)).to_int());
	//result.reserve(binomial(numeric(n+m-1), numeric(m-1)).to_int());
	intvector k(m-1);
	intvector k_cum(m-1); // k_cum[l]:=sum(i=0,l,k[l]);
	intvector upper_limit(m-1);
	int l;

	for (size_t l=0; l<m-1; ++l) {
		k[l] = 0;
		k_cum[l] = 0;
		upper_limit[l] = n;
	}

	while (true) {
		exvector term;
		term.reserve(m+1);
		for (l=0; l<m-1; ++l) {
			const ex & b = a.op(l);
			GINAC_ASSERT(!is_exactly_a<add>(b));
			GINAC_ASSERT(!is_exactly_a<power>(b) ||
			             !is_exactly_a<numeric>(ex_to<power>(b).exponent) ||
			             !ex_to<numeric>(ex_to<power>(b).exponent).is_pos_integer() ||
			             !is_exactly_a<add>(ex_to<power>(b).basis) ||
			             !is_exactly_a<mul>(ex_to<power>(b).basis) ||
			             !is_exactly_a<power>(ex_to<power>(b).basis));
			if (is_exactly_a<mul>(b))
				term.push_back(expand_mul(ex_to<mul>(b), numeric(k[l]), options, true));
			else
				term.push_back(power(b,k[l]));
		}

		const ex & b = a.op(l);
		GINAC_ASSERT(!is_exactly_a<add>(b));
		GINAC_ASSERT(!is_exactly_a<power>(b) ||
		             !is_exactly_a<numeric>(ex_to<power>(b).exponent) ||
		             !ex_to<numeric>(ex_to<power>(b).exponent).is_pos_integer() ||
		             !is_exactly_a<add>(ex_to<power>(b).basis) ||
		             !is_exactly_a<mul>(ex_to<power>(b).basis) ||
		             !is_exactly_a<power>(ex_to<power>(b).basis));
		if (is_exactly_a<mul>(b))
			term.push_back(expand_mul(ex_to<mul>(b), numeric(n-k_cum[m-2]), options, true));
		else
			term.push_back(power(b,n-k_cum[m-2]));


		numeric f = py_binomial_int(n,k[0]);
		for (l=1; l<m-1; ++l)
		  f *= py_binomial_int(n-k_cum[l-1], k[l]);

		term.push_back(f);

		result.push_back(ex((new mul(term))->setflag(status_flags::dynallocated)).expand(options));

		// increment k[]
		l = m-2;
		while ((l>=0) && ((++k[l])>upper_limit[l])) {
			k[l] = 0;
			--l;
		}
		if (l<0) break;

		// recalc k_cum[] and upper_limit[]
		k_cum[l] = (l==0 ? k[0] : k_cum[l-1]+k[l]);

		for (size_t i=l+1; i<m-1; ++i)
			k_cum[i] = k_cum[i-1]+k[i];

		for (size_t i=l+1; i<m-1; ++i)
			upper_limit[i] = n-k_cum[i-1];
	}

	return (new add(result))->setflag(status_flags::dynallocated |
	                                  status_flags::expanded);
}


/** Special case of power::expand_add. Expands a^2 where a is an add.
 *  @see power::expand_add */
ex power::expand_add_2(const add & a, unsigned options) const
{
	epvector sum;
	size_t a_nops = a.nops();
	sum.reserve((a_nops*(a_nops+1))/2);
	epvector::const_iterator last = a.seq.end();

	// power(+(x,...,z;c),2)=power(+(x,...,z;0),2)+2*c*+(x,...,z;0)+c*c
	// first part: ignore overall_coeff and expand other terms
	for (epvector::const_iterator cit0=a.seq.begin(); cit0!=last; ++cit0) {
		const ex & r = cit0->rest;
		const ex & c = cit0->coeff;
		
		GINAC_ASSERT(!is_exactly_a<add>(r));
		GINAC_ASSERT(!is_exactly_a<power>(r) ||
		             !is_exactly_a<numeric>(ex_to<power>(r).exponent) ||
		             !ex_to<numeric>(ex_to<power>(r).exponent).is_pos_integer() ||
		             !is_exactly_a<add>(ex_to<power>(r).basis) ||
		             !is_exactly_a<mul>(ex_to<power>(r).basis) ||
		             !is_exactly_a<power>(ex_to<power>(r).basis));
		
		if (c.is_equal(_ex1)) {
			if (is_exactly_a<mul>(r)) {
				sum.push_back(expair(expand_mul(ex_to<mul>(r), *_num2_p, options, true),
				                     _ex1));
			} else {
				sum.push_back(expair((new power(r,_ex2))->setflag(status_flags::dynallocated),
				                     _ex1));
			}
		} else {
			if (is_exactly_a<mul>(r)) {
				sum.push_back(a.combine_ex_with_coeff_to_pair(expand_mul(ex_to<mul>(r), *_num2_p, options, true),
				                     ex_to<numeric>(c).power_dyn(*_num2_p)));
			} else {
				sum.push_back(a.combine_ex_with_coeff_to_pair((new power(r,_ex2))->setflag(status_flags::dynallocated),
				                     ex_to<numeric>(c).power_dyn(*_num2_p)));
			}
		}

		for (epvector::const_iterator cit1=cit0+1; cit1!=last; ++cit1) {
			const ex & r1 = cit1->rest;
			const ex & c1 = cit1->coeff;
			sum.push_back(a.combine_ex_with_coeff_to_pair((new mul(r,r1))->setflag(status_flags::dynallocated),
			                                              _num2_p->mul(ex_to<numeric>(c)).mul_dyn(ex_to<numeric>(c1))));
		}
	}
	
	GINAC_ASSERT(sum.size()==(a.seq.size()*(a.seq.size()+1))/2);
	
	// second part: add terms coming from overall_factor (if != 0)
	if (!a.overall_coeff.is_zero()) {
		epvector::const_iterator i = a.seq.begin(), end = a.seq.end();
		while (i != end) {
			sum.push_back(a.combine_pair_with_coeff_to_pair(*i, ex_to<numeric>(a.overall_coeff).mul_dyn(*_num2_p)));
			++i;
		}
		sum.push_back(expair(ex_to<numeric>(a.overall_coeff).power_dyn(*_num2_p),_ex1));
	}
	
	GINAC_ASSERT(sum.size()==(a_nops*(a_nops+1))/2);
	
	return (new add(sum))->setflag(status_flags::dynallocated | status_flags::expanded);
}

/** Expand factors of m in m^n where m is a mul and n is an integer.
 *  @see power::expand */
ex power::expand_mul(const mul & m, const numeric & n, unsigned options, bool from_expand) const
{
	GINAC_ASSERT(n.is_integer());

	if (n.is_zero()) {
		return _ex1;
	}

	// do not bother to rename indices if there are no any.
	if ((!(options & expand_options::expand_rename_idx)) 
			&& m.info(info_flags::has_indices))
		options |= expand_options::expand_rename_idx;
	// Leave it to multiplication since dummy indices have to be renamed
	if ((options & expand_options::expand_rename_idx) &&
		(get_all_dummy_indices(m).size() > 0) && n.is_positive()) {
		ex result = m;
		exvector va = get_all_dummy_indices(m);
		sort(va.begin(), va.end(), ex_is_less());

		for (int i=1; i < n.to_int(); i++)
			result *= rename_dummy_indices_uniquely(va, m);
		return result;
	}

	epvector distrseq;
	distrseq.reserve(m.seq.size());
	bool need_reexpand = false;

	epvector::const_iterator last = m.seq.end();
	epvector::const_iterator cit = m.seq.begin();
	while (cit!=last) {
		expair p = m.combine_pair_with_coeff_to_pair(*cit, n);
		if (from_expand && is_exactly_a<add>(cit->rest) && ex_to<numeric>(p.coeff).is_pos_integer()) {
			// this happens when e.g. (a+b)^(1/2) gets squared and
			// the resulting product needs to be reexpanded
			need_reexpand = true;
		}
		distrseq.push_back(p);
		++cit;
	}

	const mul & result = static_cast<const mul &>((new mul(distrseq, ex_to<numeric>(m.overall_coeff).power_dyn(n)))->setflag(status_flags::dynallocated));
	if (need_reexpand)
		return ex(result).expand(options);
	if (from_expand)
		return result.setflag(status_flags::expanded);
	return result;
}

} // namespace GiNaC
