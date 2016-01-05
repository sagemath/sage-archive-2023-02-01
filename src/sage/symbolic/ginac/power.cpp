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

#include "py_funcs.h"
#include "power.h"
#include "expairseq.h"
#include "add.h"
#include "mul.h"
#include "ncmul.h"
#include "numeric.h"
#include "constant.h"
#include "infinity.h"
#include "operators.h"
#include "inifcns.h" // for log() in power::derivative() and exp for printing
#include "matrix.h"
#include "indexed.h"
#include "symbol.h"
#include "lst.h"
#include "archive.h"
#include "utils.h"
#include "relational.h"
#include "compiler.h"
#include "function.h"

#include <vector>
#include <stdexcept>
#include <limits>
#include <sstream>
#include <string>
#ifdef DO_GINAC_ASSERT
#  include <typeinfo>
#endif

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

long power::calchash() const
{
  hashvalue = basis.gethash() + exponent.gethash();
  setflag(status_flags::hash_calculated);
  return hashvalue;
}

void power::print_power(const print_context & c, const char *powersymbol, const char *openbrace, const char *closebrace, unsigned level) const
{
	// Ordinary output of powers using '^' or '**'
	if (precedence() <= level) 
		c.s << openbrace << '(';
	basis.print(c, precedence());
	if (precedence() <= level)
		c.s << ')' << closebrace;
	c.s << powersymbol;
	c.s << openbrace;
	exponent.print(c, level);
	c.s << closebrace;
}

void power::do_print_dflt(const print_dflt & c, unsigned level) const
{
       if ((-exponent).is_integer_one()) {
		// inverses printed in a special way
 	        if (level >= 20) c.s << "(";
		c.s << "1/";
		basis.print(c, precedence());
 	        if (level >= 20) c.s << ")";
	} else if (exponent.is_equal(_ex1_2)) {

		// Square roots are printed in a special way
		c.s << "sqrt(";
		basis.print(c);
		c.s << ')';

	} else if (exponent.is_equal(_ex_1_2)) {
		// Square roots are printed in a special way
		c.s << "1/sqrt(";
		basis.print(c);
		c.s << ')';
	} else {

		std::stringstream tstream;
		print_dflt tcontext(tstream, c.options);
		exponent.print(tcontext, precedence());
		std::string expstr = tstream.str();
                if (expstr[0] == '-') {
			c.s<<"1/";
			expstr = expstr.erase(0, 1);
		}
		if (precedence() <= level)
			c.s << '(';
		// exp function prints as e^a. Printing powers of this can be
		// confusing, so we add parenthesis if the basis is exp
		bool exp_parenthesis = is_ex_the_function(basis, exp) and
                        not basis.op(0).is_integer_one();
		if (exp_parenthesis)
			c.s << '(';
		basis.print(c, precedence());
		if (exp_parenthesis)
			c.s << ')';
		if (not (-exponent).is_integer_one()) {
    		        c.s << "^" << expstr;
		}
		if (precedence() <= level)
			c.s << ')';
	}
		
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

	} else {
		// exp function prints as e^a. Printing powers of this can be
		// confusing, so we add parenthesis if the basis is exp
		bool base_parenthesis = is_ex_the_function(basis, exp) and
			not basis.op(0).is_integer_one();

		if (precedence() <= level)
			c.s << "{\\left(";
		if (base_parenthesis)
			c.s << "\\left(";

		basis.print(c, precedence());

		if (base_parenthesis)
			c.s << "\\right)";

		if (not (-exponent).is_integer_one()) {
    		        c.s << "^{";
			bool exp_parenthesis = is_exactly_a<power>(exponent);
			if (exp_parenthesis)
				c.s << "\\left(";
			exponent.print(c, 0);
			if (exp_parenthesis)
				c.s << "\\right)";
			c.s << '}';
		}

		if (precedence() <= level)
			c.s << "\\right)"<<'}';
	}
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
		if (exponent.is_integer_one()) {
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
	 && (is_exactly_a<symbol>(basis) || is_exactly_a<constant>(basis))) {
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
	} else if (exponent.is_integer_one()) {
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
                case info_flags::inexact:
                        return exponent.info(info_flags::inexact) or
                                basis.info(info_flags::inexact);
		case info_flags::algebraic:
			return !exponent.info(info_flags::integer) ||
			       basis.info(inf);
		case info_flags::expanded:
			return (flags & status_flags::expanded);
		case info_flags::positive:
                        if (exponent.info(info_flags::even))
                                return basis.info(info_flags::real);
                        if (exponent.info(info_flags::odd))
                                return basis.info(info_flags::positive);
			return basis.info(info_flags::positive) && exponent.info(info_flags::real);
		case info_flags::nonnegative:
			return (basis.info(info_flags::positive) && exponent.info(info_flags::real)) ||
			       (basis.info(info_flags::real) && exponent.info(info_flags::integer) && exponent.info(info_flags::even));
		case info_flags::real:
			return basis.info(inf) && exponent.info(info_flags::integer);
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
	if (basis.is_polynomial(var)) {
		if (basis.has(var))
			// basis is non-constant polynomial in var
			return exponent.info(info_flags::nonnegint);
		else
			// basis is constant in var
			return !exponent.has(var);
	}
	// basis is a non-polynomial function of var
	return false;
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
 *  - ^(^(x,c1),c2) -> ^(x,c1*c2)  (c2 integer or -1 < c1 <= 1 or (c1=-1 and c2>0), case c1=1 should not happen, see below!)
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
	numeric num_basis;
	numeric num_exponent;
	
	if (is_exactly_a<numeric>(ebasis)) {
		basis_is_numerical = true;
		num_basis = ex_to<numeric>(ebasis);
	}
	if (is_exactly_a<numeric>(eexponent)) {
		exponent_is_numerical = true;
		num_exponent = ex_to<numeric>(eexponent);
	}

	// ^(\infty, x)
	if (is_exactly_a<infinity>(ebasis)) {
		const infinity & basis_inf = ex_to<infinity>(ebasis);
		if (eexponent.nsymbols()>0)
			throw(std::domain_error("power::eval(): pow(Infinity, f(x)) is not defined."));
		if (eexponent.is_zero())
			throw(std::domain_error("power::eval(): pow(Infinity, 0) is undefined."));
		if (eexponent.info(info_flags::negative))
			return _ex0;
		if (eexponent.info(info_flags::positive)) {
			if (basis_inf.is_unsigned_infinity())
				return UnsignedInfinity;
			else
				return mul(pow(basis_inf.get_direction(), eexponent), Infinity);
                }
		throw(std::domain_error("power::eval(): pow(Infinity, c)"
					" for constant of undetermined sign is not defined."));
	}

	// ^(x, \infty)
	if (is_exactly_a<infinity>(eexponent)) {
		const infinity & exp_inf = ex_to<infinity>(eexponent);
		if (exp_inf.is_unsigned_infinity())
			throw(std::domain_error("power::eval(): pow(x, unsigned_infinity) is not defined."));
		if (ebasis.nsymbols()>0) 
			throw(std::domain_error("power::eval(): pow(f(x), infinity) is not defined."));
		// x^(c*oo) --> (x^c)^(+oo)
		const ex abs_base = abs(pow(ebasis, exp_inf.get_direction()));
		if (abs_base > _ex1) {
			if (ebasis.info(info_flags::positive))
				return Infinity;
			else
				return UnsignedInfinity;
                        }
		if (abs_base < _ex1) return _ex0;
		if (abs_base == _ex1)
			throw(std::domain_error("power::eval(): pow(1, Infinity) is not defined."));
		throw(std::domain_error("power::eval(): pow(c, Infinity)"
					" for unknown magnitude |c| is not defined."));
	}
	
	// ^(x,0) -> 1  (0^0 also handled here)
	if (eexponent.is_zero() && 
		!(basis_is_numerical && 
			num_basis.is_parent_pos_char())) {
		return _ex1;
	}
	
	// ^(x,1) -> x
	if (eexponent.is_integer_one())
		return ebasis;

	// ^(0,c1) -> 0 or exception  (depending on real value of c1)
	if (ebasis.is_zero() && exponent_is_numerical) {
		if ((num_exponent.real()).is_zero())
			throw (std::domain_error("power::eval(): pow(0,I) is undefined"));
		else if ((num_exponent.real()).is_negative())
			throw (pole_error("power::eval(): division by zero",1));
		else
			return _ex0;
	}

	// ^(1,x) -> 1
	if (ebasis.is_integer_one())
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
			const bool basis_is_crational = num_basis.is_crational();
                        const bool exponent_is_rational = num_exponent.is_rational();
			const bool exponent_is_crational = exponent_is_rational || num_exponent.is_crational();
			if (!basis_is_crational || !exponent_is_crational) {
				// return a plain float
				return (new numeric(num_basis.power(num_exponent)))->setflag(status_flags::dynallocated |
				                                                               status_flags::evaluated |
				                                                               status_flags::expanded);
			}

			if (exponent_is_rational) {
				const numeric res = num_basis.power(num_exponent);
				if (res.is_crational()) {
					return res;
				}
			}
			GINAC_ASSERT(!num_exponent->is_integer());  // has been handled by now

			// ^(c1,n/m) -> *(c1^q,c1^(n/m-q)), 0<(n/m-q)<1, q integer
			if (basis_is_crational && exponent_is_rational) {
				const numeric n = num_exponent.numer();
				const numeric m = num_exponent.denom();
				numeric r;
				numeric q = iquo(n, m, r);
				if (r.is_negative()) {
					r += m;
					--q;
				}
				if (q.is_zero()) {  // the exponent was in the allowed range 0<(n/m)<1
					if (num_basis.is_rational()) {
						// call rational_power_parts
						// for a^b return c,d such that a^b = c*d^b
						PyObject* basis_py = num_basis.to_pyobject();
						PyObject* exponent_py = num_exponent.to_pyobject();
						PyObject* restuple = py_funcs.py_rational_power_parts(basis_py, exponent_py);
						if(!restuple) {
							throw(std::runtime_error("power::eval, error in rational_power_parts"));
						}
						PyObject *ppower = PyTuple_GetItem(restuple,0);
						PyObject *newbasis = PyTuple_GetItem(restuple,1);
						const bool ppower_equals_one = PyObject_IsTrue(PyTuple_GetItem(restuple, 2));
						ex result = (new power(newbasis,exponent))->setflag(status_flags::dynallocated | status_flags::evaluated);
						if (not ppower_equals_one)
							result = (new mul(result, ppower))->setflag(status_flags::dynallocated | status_flags::evaluated);
						Py_DECREF(restuple);
						Py_DECREF(basis_py);
						Py_DECREF(exponent_py);
						return result;
					}
					return this->hold();
				} else if (r.is_zero()) {
					// if r == 0, the following else clause causes the power
					// constructor to be called again with the same parameter
					// leading to an infinite loop
					return num_basis.power(q);
				} else {
					// assemble resulting product, but allowing for a re-evaluation,
					// because otherwise we'll end up with something like
					//    (7/8)^(4/3)  ->  7/8*(1/2*7^(1/3))
					// instead of 7/16*7^(1/3).
					ex prod = power(num_basis, r.div(m));
					return prod * power(num_basis, q);
				}
			}
		}
	
		// ^(^(x,c1),c2) -> ^(x,c1*c2)
		// (c1, c2 numeric(), c2 integer or -1 < c1 <= 1 or (c1=-1 and c2>0),
		// case c1==1 should not happen, see below!)
		if (is_exactly_a<power>(ebasis)) {
			const power & sub_power = ex_to<power>(ebasis);
			const ex & sub_basis = sub_power.basis;
			const ex & sub_exponent = sub_power.exponent;
			if (is_exactly_a<numeric>(sub_exponent)) {
				const numeric & num_sub_exponent = ex_to<numeric>(sub_exponent);
				GINAC_ASSERT(num_sub_exponent!=numeric(1));
				if (num_exponent.is_integer() || (abs(num_sub_exponent) - (*_num1_p)).is_negative() 
						|| (num_sub_exponent == *_num_1_p && num_exponent.is_positive())) {
					return power(sub_basis,num_sub_exponent.mul(num_exponent));
				}
			}
		}
	
		if (num_exponent.is_integer()) {
                        
                        // ^(*(x,y,z),c1) -> *(x^c1,y^c1,z^c1) (c1 integer)
                        if (is_exactly_a<mul>(ebasis)) {
                                return expand_mul(ex_to<mul>(ebasis), num_exponent, 0);
                        }

                        // (2*x + 6*y)^(-4) -> 1/16*(x + 3*y)^(-4)
                        if (is_exactly_a<add>(ebasis)) {
                                numeric icont = ebasis.integer_content();
                                const numeric lead_coeff = 
                                        ex_to<numeric>(ex_to<add>(ebasis).\
                                                        lead_coeff()).div(icont);

                                const bool canonicalizable = lead_coeff.is_integer();
                                const bool unit_normal = lead_coeff.is_pos_integer();
                                if (canonicalizable && (! unit_normal))
                                        icont = icont.mul(*_num_1_p);

                                if (canonicalizable && (icont != *_num1_p)) {
                                        const add& addref = ex_to<add>(ebasis);
                                        auto  addp = new add(addref);
                                        addp->setflag(status_flags::dynallocated);
                                        addp->clearflag(status_flags::hash_calculated);
                                        addp->overall_coeff = ex_to<numeric>(addp->overall_coeff).div_dyn(icont);
                                        addp->seq_sorted.resize(0);
                                        for (auto & elem : addp->seq)
                                                elem.coeff = ex_to<numeric>(elem.coeff).div_dyn(icont);

                                        const numeric c = icont.power(num_exponent);
                                        if (likely(c != *_num1_p))
                                                return (new mul(power(*addp, num_exponent), c))->setflag(status_flags::dynallocated);
                                        else
                                                return power(*addp, num_exponent);
                                }
                        }
                }

		// ^(*(...,x;c1),c2) -> *(^(*(...,x;1),c2),c1^c2)  (c1, c2 numeric(), c1>0)
		// ^(*(...,x;c1),c2) -> *(^(*(...,x;-1),c2),(-c1)^c2)  (c1, c2 numeric(), c1<0)
		if (is_exactly_a<mul>(ebasis)) {
			GINAC_ASSERT(!num_exponent.is_integer()); // should have been handled above
			const mul & mulref = ex_to<mul>(ebasis);
			if (!mulref.overall_coeff.is_equal(_ex1)) {
				const numeric & num_coeff = ex_to<numeric>(mulref.overall_coeff);
				if (num_coeff.is_real()) {
					if (num_coeff.is_positive()) {
						auto mulp = new mul(mulref);
						mulp->overall_coeff = _ex1;
						mulp->setflag(status_flags::dynallocated);
						mulp->clearflag(status_flags::evaluated);
						mulp->clearflag(status_flags::hash_calculated);
						mulp->seq_sorted.resize(0);
						return (new mul(power(*mulp,exponent),
						                power(num_coeff, num_exponent)))->setflag(status_flags::dynallocated);
					} else {
						GINAC_ASSERT(num_coeff.compare(*_num0_p)<0);
						if (!num_coeff.is_equal(*_num_1_p)) {
							auto mulp = new mul(mulref);
							mulp->overall_coeff = _ex_1;
							mulp->setflag(status_flags::dynallocated);
							mulp->clearflag(status_flags::evaluated);
							mulp->clearflag(status_flags::hash_calculated);
							mulp->seq_sorted.resize(0);
							return (new mul(power(*mulp, exponent),
							                power(abs(num_coeff), num_exponent)))->setflag(status_flags::dynallocated);
						}
					}
				}
			}
		}

		// ^(nc,c1) -> ncmul(nc,nc,...) (c1 positive integer, unless nc is a matrix)
		if (ebasis.return_type() != return_types::commutative &&
                    num_exponent.is_pos_integer() &&
                    !is_exactly_a<matrix>(ebasis)) {
			return ncmul(exvector(num_exponent.to_int(), ebasis), true);
		}
	}
	
	if (are_ex_trivially_equal(ebasis,basis) &&
	    are_ex_trivially_equal(eexponent,exponent)) {
		return this->hold();
	}
	return (new power(ebasis, eexponent))->setflag(status_flags::dynallocated |
	                                               status_flags::evaluated);
}

ex power::evalf(int level, PyObject* parent) const
{
	ex ebasis;
	ex eexponent;
	
	if (level==1) {
		ebasis = basis;
		eexponent = exponent;
	} else if (level == -max_recursion_level) {
		throw(std::runtime_error("max recursion level reached"));
	} else {
		ebasis = basis.evalf(level-1, parent);
		if (!is_exactly_a<numeric>(exponent))
			eexponent = exponent.evalf(level-1, parent);
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
			&& ex_to<numeric>(exponent)	> ex_to<numeric>(other.op(1))
			&& basis.match(other.op(0)))
		return true;
	if (exponent.info(info_flags::negint)
			&& other.op(1).info(info_flags::negint)
			&& ex_to<numeric>(exponent) < ex_to<numeric>(other.op(1))
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

        for (const auto & elem : m) {
		int nummatches = std::numeric_limits<int>::max();
		lst repls;
		if (tryfactsubs(*this, elem.first, nummatches, repls))
			return (ex_to<basic>((*this) * power(elem.second.subs(ex(repls),
                                subs_options::no_pattern) / elem.first.subs(ex(repls),
                                subs_options::no_pattern),
                                nummatches))).subs_one_level(m, options);
	}

	return subs_one_level(m, options);
}

ex power::eval_ncmul(const exvector & v) const
{
	return inherited::eval_ncmul(v);
}

ex power::conjugate() const
{
	// conjugate(pow(x,y))==pow(conjugate(x),conjugate(y)) unless on the
	// branch cut which runs along the negative real axis.
	if (basis.info(info_flags::positive)) {
		ex newexponent = exponent.conjugate();
		if (are_ex_trivially_equal(exponent, newexponent)) {
			return *this;
		}
		return (new power(basis, newexponent))->setflag(status_flags::dynallocated);
	}
	if (exponent.info(info_flags::integer)) {
		ex newbasis = basis.conjugate();
		if (are_ex_trivially_equal(basis, newbasis)) {
			return *this;
		}
		return (new power(newbasis, exponent))->setflag(status_flags::dynallocated);
	}
	return conjugate_function(*this).hold();
}

ex power::real_part() const
{
	// basis == a+I*b, exponent == c+I*d
	const ex a = basis.real_part();
	const ex c = exponent.real_part();
	if (basis.is_equal(a) && exponent.is_equal(c)) {
		// Re(a^c)
		return *this;
	}

	const ex b = basis.imag_part();
	if (exponent.info(info_flags::integer)) {
		// Re((a+I*b)^c)  w/  c ∈ ℤ
		long N = ex_to<numeric>(c).to_long();
		// Use real terms in Binomial expansion to construct
		// Re(expand(power(a+I*b, N))).
		long NN = N > 0 ? N : -N;
		ex numer = N > 0 ? _ex1 : power(power(a,2) + power(b,2), NN);
		ex result = 0;
		for (long n = 0; n <= NN; n += 2) {
			ex term = binomial(NN, n) * power(a, NN-n) * power(b, n) / numer;
			if (n % 4 == 0) {
				result += term;  // sign: I^n w/ n == 4*m
			} else {
				result -= term;  // sign: I^n w/ n == 4*m+2
			}
		}
		return result;
	}

	// Re((a+I*b)^(c+I*d))
	const ex d = exponent.imag_part();
	return power(abs(basis),c)*exp(-d*atan2(b,a))*cos(c*atan2(b,a)+d*log(abs(basis)));
}

ex power::imag_part() const
{
	const ex a = basis.real_part();
	const ex c = exponent.real_part();
	if (basis.is_equal(a) && exponent.is_equal(c)) {
		// Im(a^c)
		return 0;
	}

	const ex b = basis.imag_part();
	if (exponent.info(info_flags::integer)) {
		// Im((a+I*b)^c)  w/  c ∈ ℤ
		long N = ex_to<numeric>(c).to_long();
		// Use imaginary terms in Binomial expansion to construct
		// Im(expand(power(a+I*b, N))).
		long p = N > 0 ? 1 : 3;  // modulus for positive sign
		long NN = N > 0 ? N : -N;
		ex numer = N > 0 ? _ex1 : power(power(a,2) + power(b,2), NN);
		ex result = 0;
		for (long n = 1; n <= NN; n += 2) {
			ex term = binomial(NN, n) * power(a, NN-n) * power(b, n) / numer;
			if (n % 4 == p) {
				result += term;  // sign: I^n w/ n == 4*m+p
			} else {
				result -= term;  // sign: I^n w/ n == 4*m+2+p
			}
		}
		return result;
	}

	// Im((a+I*b)^(c+I*d))
	const ex d = exponent.imag_part();
	return power(abs(basis),c)*exp(-d*atan2(b,a))*sin(c*atan2(b,a)+d*log(abs(basis)));
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
	    // If the exponent is not a function of s, we have the following nice
	    // looking formula.  We use this to avoid getting ugly and hard to
	    // read output.
	    // D(b^e) = e * b^(e-1) * D(b)
	    ex ediff = exponent.diff(s);
	    if (ediff == 0) 
		return mul(mul(exponent, power(basis, exponent-1)), basis.diff(s));
	    
	    // The formula in the general case.
	    // D(b^e) = b^e * (D(e)*ln(b) + e*D(b)/b)
	    return mul(*this,
		       add(mul(ediff, log(basis)),
			   mul(mul(exponent, basis.diff(s)), power(basis, _ex_1))));
	}
}

int power::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<power>(other));
	const power &o = static_cast<const power &>(other);

	int cmpval = basis.compare(o.basis);
	if (cmpval)
		return cmpval;
	else
		return exponent.compare(o.exponent);
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

	// (x*p)^c -> x^c * p^c, if p>0
	// makes sense before expanding the basis
	if (is_exactly_a<mul>(basis) && !basis.info(info_flags::indefinite)) {
		const mul &m = ex_to<mul>(basis);
		exvector prodseq;
		epvector powseq;
		prodseq.reserve(m.seq.size() + 1);
		powseq.reserve(m.seq.size() + 1);
		bool possign = true;

		// search for positive/negative factors
                for (const auto & elem : m.seq) {
			const ex& e = m.recombine_pair_to_ex(elem);
			if (e.info(info_flags::positive))
				prodseq.push_back(pow(e, exponent).expand(options));
			else if (e.info(info_flags::negative)) {
				prodseq.push_back(pow(-e, exponent).expand(options));
				possign = !possign;
			} else
				powseq.push_back(elem);
		}

		// take care on the numeric coefficient
		ex coef = (possign? _ex1 : _ex_1);
		if (m.overall_coeff.info(info_flags::positive) && m.overall_coeff != _ex1)
			prodseq.push_back(power(m.overall_coeff, exponent));
		else if (m.overall_coeff.info(info_flags::negative) && m.overall_coeff != _ex_1)
			prodseq.push_back(power(-m.overall_coeff, exponent));
		else
			coef *= m.overall_coeff;

		// If positive/negative factors are found, then extract them.
		// In either case we set a flag to avoid the second run on a part
		// which does not have positive/negative terms.
		if (prodseq.size() > 0) {
			ex newbasis = coef*mul(powseq);
			ex_to<basic>(newbasis).setflag(status_flags::purely_indefinite);
			return ((new mul(prodseq))->setflag(status_flags::dynallocated)*(new power(newbasis, exponent))->setflag(status_flags::dynallocated).expand(options)).expand(options);
		} else
			ex_to<basic>(basis).setflag(status_flags::purely_indefinite);
	}

	const ex expanded_basis = basis.expand(options);
	const ex expanded_exponent = exponent.expand(options);
	
	// x^(a+b) -> x^a * x^b
	if (is_exactly_a<add>(expanded_exponent)) {
		const add &a = ex_to<add>(expanded_exponent);
		exvector distrseq;
		distrseq.reserve(a.seq.size() + 1);
                for (const auto & elem : a.seq)
			distrseq.push_back(power(expanded_basis, a.recombine_pair_to_ex(elem)));
		
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
	result.reserve(numeric(py_funcs.py_binomial_int(n+m-1, m-1)).to_int());
	//result.reserve(binomial(numeric(n+m-1), numeric(m-1)).to_int());
	intvector k(m-1);
	intvector k_cum(m-1); // k_cum[l]:=sum(i=0,l,k[l]);
	intvector upper_limit(m-1);

	for (size_t l=0; l<m-1; ++l) {
		k[l] = 0;
		k_cum[l] = 0;
		upper_limit[l] = n;
	}

	while (true) {
		exvector term;
		term.reserve(m+1);
		for (size_t l=0; l<m-1; ++l) {
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

		const ex & b = a.op(m-1);
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


		numeric f = py_funcs.py_binomial_int(n,k[0]);
		for (size_t l=1; l<m-1; ++l)
		  f *= py_funcs.py_binomial_int(n-k_cum[l-1], k[l]);

		term.push_back(f);

		result.push_back(ex((new mul(term))->setflag(status_flags::dynallocated)).expand(options));

		// increment k[]
		int l = m-2;
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
	auto last = a.seq.end();

	// power(+(x,...,z;c),2)=power(+(x,...,z;0),2)+2*c*+(x,...,z;0)+c*c
	// first part: ignore overall_coeff and expand other terms
	for (auto cit0=a.seq.begin(); cit0!=last; ++cit0) {
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
				sum.push_back(a.combine_ex_with_coeff_to_pair(expand_mul(ex_to<mul>(r), *_num2_p, options, true),
				                                              _ex1));
			} else {
				sum.push_back(a.combine_ex_with_coeff_to_pair((new power(r,_ex2))->setflag(status_flags::dynallocated),
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

		for (auto cit1=cit0+1; cit1!=last; ++cit1) {
			const ex & r1 = cit1->rest;
			const ex & c1 = cit1->coeff;
			sum.push_back(a.combine_ex_with_coeff_to_pair(mul(r,r1).expand(options),
			                                              _num2_p->mul(ex_to<numeric>(c)).mul_dyn(ex_to<numeric>(c1))));
		}
	}
	
	GINAC_ASSERT(sum.size()==(a.seq.size()*(a.seq.size()+1))/2);
	
	// second part: add terms coming from overall_coeff (if != 0)
	if (!a.overall_coeff.is_zero()) {
                for (const auto & elem : a.seq)
			sum.push_back(a.combine_pair_with_coeff_to_pair(elem, ex_to<numeric>(a.overall_coeff).mul_dyn(*_num2_p)));
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

        for (const auto & elem : m.seq) {
		expair p = m.combine_pair_with_coeff_to_pair(elem, n);
		if (from_expand && is_exactly_a<add>(elem.rest) && ex_to<numeric>(p.coeff).is_pos_integer()) {
			// this happens when e.g. (a+b)^(1/2) gets squared and
			// the resulting product needs to be reexpanded
			need_reexpand = true;
		}
		distrseq.push_back(p);
	}

	const mul & result = static_cast<const mul &>((new mul(distrseq, ex_to<numeric>(m.overall_coeff).power_dyn(n)))->setflag(status_flags::dynallocated));
	if (need_reexpand)
		return ex(result).expand(options);
	if (from_expand)
		return result.setflag(status_flags::expanded);
	return result;
}

} // namespace GiNaC
