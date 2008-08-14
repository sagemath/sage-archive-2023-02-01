/** @file normal.cpp
 *
 *  This file implements several functions that work on univariate and
 *  multivariate polynomials and rational functions.
 *  These functions include polynomial quotient and remainder, GCD and LCM
 *  computation, square-free factorization and rational function normalization. */

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
#include <map>

#include "normal.h"
#include "basic.h"
#include "ex.h"
#include "add.h"
#include "constant.h"
#include "expairseq.h"
#include "fail.h"
#include "inifcns.h"
#include "lst.h"
#include "mul.h"
#include "numeric.h"
#include "power.h"
#include "relational.h"
#include "operators.h"
#include "matrix.h"
#include "pseries.h"
#include "symbol.h"
#include "utils.h"

namespace GiNaC {

// If comparing expressions (ex::compare()) is fast, you can set this to 1.
// Some routines like quo(), rem() and gcd() will then return a quick answer
// when they are called with two identical arguments.
#define FAST_COMPARE 1

// Set this if you want divide_in_z() to use remembering
#define USE_REMEMBER 0

// Set this if you want divide_in_z() to use trial division followed by
// polynomial interpolation (always slower except for completely dense
// polynomials)
#define USE_TRIAL_DIVISION 0

// Set this to enable some statistical output for the GCD routines
#define STATISTICS 0


#if STATISTICS
// Statistics variables
static int gcd_called = 0;
static int sr_gcd_called = 0;
static int heur_gcd_called = 0;
static int heur_gcd_failed = 0;

// Print statistics at end of program
static struct _stat_print {
	_stat_print() {}
	~_stat_print() {
		std::cout << "gcd() called " << gcd_called << " times\n";
		std::cout << "sr_gcd() called " << sr_gcd_called << " times\n";
		std::cout << "heur_gcd() called " << heur_gcd_called << " times\n";
		std::cout << "heur_gcd() failed " << heur_gcd_failed << " times\n";
	}
} stat_print;
#endif


/** Return pointer to first symbol found in expression.  Due to GiNaC's
 *  internal ordering of terms, it may not be obvious which symbol this
 *  function returns for a given expression.
 *
 *  @param e  expression to search
 *  @param x  first symbol found (returned)
 *  @return "false" if no symbol was found, "true" otherwise */
static bool get_first_symbol(const ex &e, ex &x)
{
	if (is_a<symbol>(e)) {
		x = e;
		return true;
	} else if (is_exactly_a<add>(e) || is_exactly_a<mul>(e)) {
		for (size_t i=0; i<e.nops(); i++)
			if (get_first_symbol(e.op(i), x))
				return true;
	} else if (is_exactly_a<power>(e)) {
		if (get_first_symbol(e.op(0), x))
			return true;
	}
	return false;
}


/*
 *  Statistical information about symbols in polynomials
 */

/** This structure holds information about the highest and lowest degrees
 *  in which a symbol appears in two multivariate polynomials "a" and "b".
 *  A vector of these structures with information about all symbols in
 *  two polynomials can be created with the function get_symbol_stats().
 *
 *  @see get_symbol_stats */
struct sym_desc {
	/** Reference to symbol */
	ex sym;

	/** Highest degree of symbol in polynomial "a" */
	int deg_a;

	/** Highest degree of symbol in polynomial "b" */
	int deg_b;

	/** Lowest degree of symbol in polynomial "a" */
	int ldeg_a;

	/** Lowest degree of symbol in polynomial "b" */
	int ldeg_b;

	/** Maximum of deg_a and deg_b (Used for sorting) */
	int max_deg;

	/** Maximum number of terms of leading coefficient of symbol in both polynomials */
	size_t max_lcnops;

	/** Commparison operator for sorting */
	bool operator<(const sym_desc &x) const
	{
		if (max_deg == x.max_deg)
			return max_lcnops < x.max_lcnops;
		else
			return max_deg < x.max_deg;
	}
};

// Vector of sym_desc structures
typedef std::vector<sym_desc> sym_desc_vec;

// Add symbol the sym_desc_vec (used internally by get_symbol_stats())
static void add_symbol(const ex &s, sym_desc_vec &v)
{
	sym_desc_vec::const_iterator it = v.begin(), itend = v.end();
	while (it != itend) {
		if (it->sym.is_equal(s))  // If it's already in there, don't add it a second time
			return;
		++it;
	}
	sym_desc d;
	d.sym = s;
	v.push_back(d);
}

// Collect all symbols of an expression (used internally by get_symbol_stats())
static void collect_symbols(const ex &e, sym_desc_vec &v)
{
	if (is_a<symbol>(e)) {
		add_symbol(e, v);
	} else if (is_exactly_a<add>(e) || is_exactly_a<mul>(e)) {
		for (size_t i=0; i<e.nops(); i++)
			collect_symbols(e.op(i), v);
	} else if (is_exactly_a<power>(e)) {
		collect_symbols(e.op(0), v);
	}
}

/** Collect statistical information about symbols in polynomials.
 *  This function fills in a vector of "sym_desc" structs which contain
 *  information about the highest and lowest degrees of all symbols that
 *  appear in two polynomials. The vector is then sorted by minimum
 *  degree (lowest to highest). The information gathered by this
 *  function is used by the GCD routines to identify trivial factors
 *  and to determine which variable to choose as the main variable
 *  for GCD computation.
 *
 *  @param a  first multivariate polynomial
 *  @param b  second multivariate polynomial
 *  @param v  vector of sym_desc structs (filled in) */
static void get_symbol_stats(const ex &a, const ex &b, sym_desc_vec &v)
{
	collect_symbols(a.eval(), v);   // eval() to expand assigned symbols
	collect_symbols(b.eval(), v);
	sym_desc_vec::iterator it = v.begin(), itend = v.end();
	while (it != itend) {
		int deg_a = a.degree(it->sym);
		int deg_b = b.degree(it->sym);
		it->deg_a = deg_a;
		it->deg_b = deg_b;
		it->max_deg = std::max(deg_a, deg_b);
		it->max_lcnops = std::max(a.lcoeff(it->sym).nops(), b.lcoeff(it->sym).nops());
		it->ldeg_a = a.ldegree(it->sym);
		it->ldeg_b = b.ldegree(it->sym);
		++it;
	}
	std::sort(v.begin(), v.end());

#if 0
	std::clog << "Symbols:\n";
	it = v.begin(); itend = v.end();
	while (it != itend) {
		std::clog << " " << it->sym << ": deg_a=" << it->deg_a << ", deg_b=" << it->deg_b << ", ldeg_a=" << it->ldeg_a << ", ldeg_b=" << it->ldeg_b << ", max_deg=" << it->max_deg << ", max_lcnops=" << it->max_lcnops << endl;
		std::clog << "  lcoeff_a=" << a.lcoeff(it->sym) << ", lcoeff_b=" << b.lcoeff(it->sym) << endl;
		++it;
	}
#endif
}


/*
 *  Computation of LCM of denominators of coefficients of a polynomial
 */

// Compute LCM of denominators of coefficients by going through the
// expression recursively (used internally by lcm_of_coefficients_denominators())
static numeric lcmcoeff(const ex &e, const numeric &l)
{
	if (e.info(info_flags::rational))
		return lcm(ex_to<numeric>(e).denom(), l);
	else if (is_exactly_a<add>(e)) {
		numeric c = *_num1_p;
		for (size_t i=0; i<e.nops(); i++)
			c = lcmcoeff(e.op(i), c);
		return lcm(c, l);
	} else if (is_exactly_a<mul>(e)) {
		numeric c = *_num1_p;
		for (size_t i=0; i<e.nops(); i++)
			c *= lcmcoeff(e.op(i), *_num1_p);
		return lcm(c, l);
	} else if (is_exactly_a<power>(e)) {
		if (is_a<symbol>(e.op(0)))
			return l;
		else
			return pow(lcmcoeff(e.op(0), l), ex_to<numeric>(e.op(1)));
	}
	return l;
}

/** Compute LCM of denominators of coefficients of a polynomial.
 *  Given a polynomial with rational coefficients, this function computes
 *  the LCM of the denominators of all coefficients. This can be used
 *  to bring a polynomial from Q[X] to Z[X].
 *
 *  @param e  multivariate polynomial (need not be expanded)
 *  @return LCM of denominators of coefficients */
static numeric lcm_of_coefficients_denominators(const ex &e)
{
	return lcmcoeff(e, *_num1_p);
}

/** Bring polynomial from Q[X] to Z[X] by multiplying in the previously
 *  determined LCM of the coefficient's denominators.
 *
 *  @param e  multivariate polynomial (need not be expanded)
 *  @param lcm  LCM to multiply in */
static ex multiply_lcm(const ex &e, const numeric &lcm)
{
	if (is_exactly_a<mul>(e)) {
		size_t num = e.nops();
		exvector v; v.reserve(num + 1);
		numeric lcm_accum = *_num1_p;
		for (size_t i=0; i<num; i++) {
			numeric op_lcm = lcmcoeff(e.op(i), *_num1_p);
			v.push_back(multiply_lcm(e.op(i), op_lcm));
			lcm_accum *= op_lcm;
		}
		v.push_back(lcm / lcm_accum);
		return (new mul(v))->setflag(status_flags::dynallocated);
	} else if (is_exactly_a<add>(e)) {
		size_t num = e.nops();
		exvector v; v.reserve(num);
		for (size_t i=0; i<num; i++)
			v.push_back(multiply_lcm(e.op(i), lcm));
		return (new add(v))->setflag(status_flags::dynallocated);
	} else if (is_exactly_a<power>(e)) {
		if (is_a<symbol>(e.op(0)))
			return e * lcm;
		else
			return pow(multiply_lcm(e.op(0), lcm.power(ex_to<numeric>(e.op(1)).inverse())), e.op(1));
	} else
		return e * lcm;
}


/** Compute the integer content (= GCD of all numeric coefficients) of an
 *  expanded polynomial. For a polynomial with rational coefficients, this
 *  returns g/l where g is the GCD of the coefficients' numerators and l
 *  is the LCM of the coefficients' denominators.
 *
 *  @return integer content */
numeric ex::integer_content() const
{
	return bp->integer_content();
}

numeric basic::integer_content() const
{
	return *_num1_p;
}

numeric numeric::integer_content() const
{
	return abs(*this);
}

numeric add::integer_content() const
{
	epvector::const_iterator it = seq.begin();
	epvector::const_iterator itend = seq.end();
	numeric c = *_num0_p, l = *_num1_p;
	while (it != itend) {
		GINAC_ASSERT(!is_exactly_a<numeric>(it->rest));
		GINAC_ASSERT(is_exactly_a<numeric>(it->coeff));
		c = gcd(ex_to<numeric>(it->coeff).numer(), c);
		l = lcm(ex_to<numeric>(it->coeff).denom(), l);
		it++;
	}
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	c = gcd(ex_to<numeric>(overall_coeff).numer(), c);
	l = lcm(ex_to<numeric>(overall_coeff).denom(), l);
	return c/l;
}

numeric mul::integer_content() const
{
#ifdef DO_GINAC_ASSERT
	epvector::const_iterator it = seq.begin();
	epvector::const_iterator itend = seq.end();
	while (it != itend) {
		GINAC_ASSERT(!is_exactly_a<numeric>(recombine_pair_to_ex(*it)));
		++it;
	}
#endif // def DO_GINAC_ASSERT
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	return abs(ex_to<numeric>(overall_coeff));
}


/*
 *  Polynomial quotients and remainders
 */

/** Quotient q(x) of polynomials a(x) and b(x) in Q[x].
 *  It satisfies a(x)=b(x)*q(x)+r(x).
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return quotient of a and b in Q[x] */
ex quo(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("quo: division by zero"));
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b))
		return a / b;
#if FAST_COMPARE
	if (a.is_equal(b))
		return _ex1;
#endif
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("quo: arguments must be polynomials over the rationals"));

	// Polynomial long division
	ex r = a.expand();
	if (r.is_zero())
		return r;
	int bdeg = b.degree(x);
	int rdeg = r.degree(x);
	ex blcoeff = b.expand().coeff(x, bdeg);
	bool blcoeff_is_numeric = is_exactly_a<numeric>(blcoeff);
	exvector v; v.reserve(std::max(rdeg - bdeg + 1, 0));
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (blcoeff_is_numeric)
			term = rcoeff / blcoeff;
		else {
			if (!divide(rcoeff, blcoeff, term, false))
				return (new fail())->setflag(status_flags::dynallocated);
		}
		term *= power(x, rdeg - bdeg);
		v.push_back(term);
		r -= (term * b).expand();
		if (r.is_zero())
			break;
		rdeg = r.degree(x);
	}
	return (new add(v))->setflag(status_flags::dynallocated);
}


/** Remainder r(x) of polynomials a(x) and b(x) in Q[x].
 *  It satisfies a(x)=b(x)*q(x)+r(x).
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return remainder of a(x) and b(x) in Q[x] */
ex rem(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("rem: division by zero"));
	if (is_exactly_a<numeric>(a)) {
		if  (is_exactly_a<numeric>(b))
			return _ex0;
		else
			return a;
	}
#if FAST_COMPARE
	if (a.is_equal(b))
		return _ex0;
#endif
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("rem: arguments must be polynomials over the rationals"));

	// Polynomial long division
	ex r = a.expand();
	if (r.is_zero())
		return r;
	int bdeg = b.degree(x);
	int rdeg = r.degree(x);
	ex blcoeff = b.expand().coeff(x, bdeg);
	bool blcoeff_is_numeric = is_exactly_a<numeric>(blcoeff);
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (blcoeff_is_numeric)
			term = rcoeff / blcoeff;
		else {
			if (!divide(rcoeff, blcoeff, term, false))
				return (new fail())->setflag(status_flags::dynallocated);
		}
		term *= power(x, rdeg - bdeg);
		r -= (term * b).expand();
		if (r.is_zero())
			break;
		rdeg = r.degree(x);
	}
	return r;
}


/** Decompose rational function a(x)=N(x)/D(x) into P(x)+n(x)/D(x)
 *  with degree(n, x) < degree(D, x).
 *
 *  @param a rational function in x
 *  @param x a is a function of x
 *  @return decomposed function. */
ex decomp_rational(const ex &a, const ex &x)
{
	ex nd = numer_denom(a);
	ex numer = nd.op(0), denom = nd.op(1);
	ex q = quo(numer, denom, x);
	if (is_exactly_a<fail>(q))
		return a;
	else
		return q + rem(numer, denom, x) / denom;
}


/** Pseudo-remainder of polynomials a(x) and b(x) in Q[x].
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return pseudo-remainder of a(x) and b(x) in Q[x] */
ex prem(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("prem: division by zero"));
	if (is_exactly_a<numeric>(a)) {
		if (is_exactly_a<numeric>(b))
			return _ex0;
		else
			return b;
	}
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("prem: arguments must be polynomials over the rationals"));

	// Polynomial long division
	ex r = a.expand();
	ex eb = b.expand();
	int rdeg = r.degree(x);
	int bdeg = eb.degree(x);
	ex blcoeff;
	if (bdeg <= rdeg) {
		blcoeff = eb.coeff(x, bdeg);
		if (bdeg == 0)
			eb = _ex0;
		else
			eb -= blcoeff * power(x, bdeg);
	} else
		blcoeff = _ex1;

	int delta = rdeg - bdeg + 1, i = 0;
	while (rdeg >= bdeg && !r.is_zero()) {
		ex rlcoeff = r.coeff(x, rdeg);
		ex term = (power(x, rdeg - bdeg) * eb * rlcoeff).expand();
		if (rdeg == 0)
			r = _ex0;
		else
			r -= rlcoeff * power(x, rdeg);
		r = (blcoeff * r).expand() - term;
		rdeg = r.degree(x);
		i++;
	}
	return power(blcoeff, delta - i) * r;
}


/** Sparse pseudo-remainder of polynomials a(x) and b(x) in Q[x].
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return sparse pseudo-remainder of a(x) and b(x) in Q[x] */
ex sprem(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("prem: division by zero"));
	if (is_exactly_a<numeric>(a)) {
		if (is_exactly_a<numeric>(b))
			return _ex0;
		else
			return b;
	}
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("prem: arguments must be polynomials over the rationals"));

	// Polynomial long division
	ex r = a.expand();
	ex eb = b.expand();
	int rdeg = r.degree(x);
	int bdeg = eb.degree(x);
	ex blcoeff;
	if (bdeg <= rdeg) {
		blcoeff = eb.coeff(x, bdeg);
		if (bdeg == 0)
			eb = _ex0;
		else
			eb -= blcoeff * power(x, bdeg);
	} else
		blcoeff = _ex1;

	while (rdeg >= bdeg && !r.is_zero()) {
		ex rlcoeff = r.coeff(x, rdeg);
		ex term = (power(x, rdeg - bdeg) * eb * rlcoeff).expand();
		if (rdeg == 0)
			r = _ex0;
		else
			r -= rlcoeff * power(x, rdeg);
		r = (blcoeff * r).expand() - term;
		rdeg = r.degree(x);
	}
	return r;
}


/** Exact polynomial division of a(X) by b(X) in Q[X].
 *  
 *  @param a  first multivariate polynomial (dividend)
 *  @param b  second multivariate polynomial (divisor)
 *  @param q  quotient (returned)
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return "true" when exact division succeeds (quotient returned in q),
 *          "false" otherwise (q left untouched) */
bool divide(const ex &a, const ex &b, ex &q, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("divide: division by zero"));
	if (a.is_zero()) {
		q = _ex0;
		return true;
	}
	if (is_exactly_a<numeric>(b)) {
		q = a / b;
		return true;
	} else if (is_exactly_a<numeric>(a))
		return false;
#if FAST_COMPARE
	if (a.is_equal(b)) {
		q = _ex1;
		return true;
	}
#endif
	if (check_args && (!a.info(info_flags::rational_polynomial) ||
	                   !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("divide: arguments must be polynomials over the rationals"));

	// Find first symbol
	ex x;
	if (!get_first_symbol(a, x) && !get_first_symbol(b, x))
		throw(std::invalid_argument("invalid expression in divide()"));

	// Try to avoid expanding partially factored expressions.
	if (is_exactly_a<mul>(b)) {
	// Divide sequentially by each term
		ex rem_new, rem_old = a;
		for (size_t i=0; i < b.nops(); i++) {
			if (! divide(rem_old, b.op(i), rem_new, false))
				return false;
			rem_old = rem_new;
		}
		q = rem_new;
		return true;
	} else if (is_exactly_a<power>(b)) {
		const ex& bb(b.op(0));
		int exp_b = ex_to<numeric>(b.op(1)).to_int();
		ex rem_new, rem_old = a;
		for (int i=exp_b; i>0; i--) {
			if (! divide(rem_old, bb, rem_new, false))
				return false;
			rem_old = rem_new;
		}
		q = rem_new;
		return true;
	} 
	
	if (is_exactly_a<mul>(a)) {
		// Divide sequentially each term. If some term in a is divisible 
		// by b we are done... and if not, we can't really say anything.
		size_t i;
		ex rem_i;
		bool divisible_p = false;
		for (i=0; i < a.nops(); ++i) {
			if (divide(a.op(i), b, rem_i, false)) {
				divisible_p = true;
				break;
			}
		}
		if (divisible_p) {
			exvector resv;
			resv.reserve(a.nops());
			for (size_t j=0; j < a.nops(); j++) {
				if (j==i)
					resv.push_back(rem_i);
				else
					resv.push_back(a.op(j));
			}
			q = (new mul(resv))->setflag(status_flags::dynallocated);
			return true;
		}
	} else if (is_exactly_a<power>(a)) {
		// The base itself might be divisible by b, in that case we don't
		// need to expand a
		const ex& ab(a.op(0));
		int a_exp = ex_to<numeric>(a.op(1)).to_int();
		ex rem_i;
		if (divide(ab, b, rem_i, false)) {
			q = rem_i*power(ab, a_exp - 1);
			return true;
		}
		for (int i=2; i < a_exp; i++) {
			if (divide(power(ab, i), b, rem_i, false)) {
				q = rem_i*power(ab, a_exp - i);
				return true;
			}
		} // ... so we *really* need to expand expression.
	}
	
	// Polynomial long division (recursive)
	ex r = a.expand();
	if (r.is_zero()) {
		q = _ex0;
		return true;
	}
	int bdeg = b.degree(x);
	int rdeg = r.degree(x);
	ex blcoeff = b.expand().coeff(x, bdeg);
	bool blcoeff_is_numeric = is_exactly_a<numeric>(blcoeff);
	exvector v; v.reserve(std::max(rdeg - bdeg + 1, 0));
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (blcoeff_is_numeric)
			term = rcoeff / blcoeff;
		else
			if (!divide(rcoeff, blcoeff, term, false))
				return false;
		term *= power(x, rdeg - bdeg);
		v.push_back(term);
		r -= (term * b).expand();
		if (r.is_zero()) {
			q = (new add(v))->setflag(status_flags::dynallocated);
			return true;
		}
		rdeg = r.degree(x);
	}
	return false;
}


#if USE_REMEMBER
/*
 *  Remembering
 */

typedef std::pair<ex, ex> ex2;
typedef std::pair<ex, bool> exbool;

struct ex2_less {
	bool operator() (const ex2 &p, const ex2 &q) const 
	{
		int cmp = p.first.compare(q.first);
		return ((cmp<0) || (!(cmp>0) && p.second.compare(q.second)<0));
	}
};

typedef std::map<ex2, exbool, ex2_less> ex2_exbool_remember;
#endif


/** Exact polynomial division of a(X) by b(X) in Z[X].
 *  This functions works like divide() but the input and output polynomials are
 *  in Z[X] instead of Q[X] (i.e. they have integer coefficients). Unlike
 *  divide(), it doesn't check whether the input polynomials really are integer
 *  polynomials, so be careful of what you pass in. Also, you have to run
 *  get_symbol_stats() over the input polynomials before calling this function
 *  and pass an iterator to the first element of the sym_desc vector. This
 *  function is used internally by the heur_gcd().
 *  
 *  @param a  first multivariate polynomial (dividend)
 *  @param b  second multivariate polynomial (divisor)
 *  @param q  quotient (returned)
 *  @param var  iterator to first element of vector of sym_desc structs
 *  @return "true" when exact division succeeds (the quotient is returned in
 *          q), "false" otherwise.
 *  @see get_symbol_stats, heur_gcd */
static bool divide_in_z(const ex &a, const ex &b, ex &q, sym_desc_vec::const_iterator var)
{
	q = _ex0;
	if (b.is_zero())
		throw(std::overflow_error("divide_in_z: division by zero"));
	if (b.is_equal(_ex1)) {
		q = a;
		return true;
	}
	if (is_exactly_a<numeric>(a)) {
		if (is_exactly_a<numeric>(b)) {
			q = a / b;
			return q.info(info_flags::integer);
		} else
			return false;
	}
#if FAST_COMPARE
	if (a.is_equal(b)) {
		q = _ex1;
		return true;
	}
#endif

#if USE_REMEMBER
	// Remembering
	static ex2_exbool_remember dr_remember;
	ex2_exbool_remember::const_iterator remembered = dr_remember.find(ex2(a, b));
	if (remembered != dr_remember.end()) {
		q = remembered->second.first;
		return remembered->second.second;
	}
#endif

	if (is_exactly_a<power>(b)) {
		const ex& bb(b.op(0));
		ex qbar = a;
		int exp_b = ex_to<numeric>(b.op(1)).to_int();
		for (int i=exp_b; i>0; i--) {
			if (!divide_in_z(qbar, bb, q, var))
				return false;
			qbar = q;
		}
		return true;
	}

	if (is_exactly_a<mul>(b)) {
		ex qbar = a;
		for (const_iterator itrb = b.begin(); itrb != b.end(); ++itrb) {
			sym_desc_vec sym_stats;
			get_symbol_stats(a, *itrb, sym_stats);
			if (!divide_in_z(qbar, *itrb, q, sym_stats.begin()))
				return false;

			qbar = q;
		}
		return true;
	}

	// Main symbol
	const ex &x = var->sym;

	// Compare degrees
	int adeg = a.degree(x), bdeg = b.degree(x);
	if (bdeg > adeg)
		return false;

#if USE_TRIAL_DIVISION

	// Trial division with polynomial interpolation
	int i, k;

	// Compute values at evaluation points 0..adeg
	vector<numeric> alpha; alpha.reserve(adeg + 1);
	exvector u; u.reserve(adeg + 1);
	numeric point = *_num0_p;
	ex c;
	for (i=0; i<=adeg; i++) {
		ex bs = b.subs(x == point, subs_options::no_pattern);
		while (bs.is_zero()) {
			point += *_num1_p;
			bs = b.subs(x == point, subs_options::no_pattern);
		}
		if (!divide_in_z(a.subs(x == point, subs_options::no_pattern), bs, c, var+1))
			return false;
		alpha.push_back(point);
		u.push_back(c);
		point += *_num1_p;
	}

	// Compute inverses
	vector<numeric> rcp; rcp.reserve(adeg + 1);
	rcp.push_back(*_num0_p);
	for (k=1; k<=adeg; k++) {
		numeric product = alpha[k] - alpha[0];
		for (i=1; i<k; i++)
			product *= alpha[k] - alpha[i];
		rcp.push_back(product.inverse());
	}

	// Compute Newton coefficients
	exvector v; v.reserve(adeg + 1);
	v.push_back(u[0]);
	for (k=1; k<=adeg; k++) {
		ex temp = v[k - 1];
		for (i=k-2; i>=0; i--)
			temp = temp * (alpha[k] - alpha[i]) + v[i];
		v.push_back((u[k] - temp) * rcp[k]);
	}

	// Convert from Newton form to standard form
	c = v[adeg];
	for (k=adeg-1; k>=0; k--)
		c = c * (x - alpha[k]) + v[k];

	if (c.degree(x) == (adeg - bdeg)) {
		q = c.expand();
		return true;
	} else
		return false;

#else

	// Polynomial long division (recursive)
	ex r = a.expand();
	if (r.is_zero())
		return true;
	int rdeg = adeg;
	ex eb = b.expand();
	ex blcoeff = eb.coeff(x, bdeg);
	exvector v; v.reserve(std::max(rdeg - bdeg + 1, 0));
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (!divide_in_z(rcoeff, blcoeff, term, var+1))
			break;
		term = (term * power(x, rdeg - bdeg)).expand();
		v.push_back(term);
		r -= (term * eb).expand();
		if (r.is_zero()) {
			q = (new add(v))->setflag(status_flags::dynallocated);
#if USE_REMEMBER
			dr_remember[ex2(a, b)] = exbool(q, true);
#endif
			return true;
		}
		rdeg = r.degree(x);
	}
#if USE_REMEMBER
	dr_remember[ex2(a, b)] = exbool(q, false);
#endif
	return false;

#endif
}


/*
 *  Separation of unit part, content part and primitive part of polynomials
 */

/** Compute unit part (= sign of leading coefficient) of a multivariate
 *  polynomial in Q[x]. The product of unit part, content part, and primitive
 *  part is the polynomial itself.
 *
 *  @param x  main variable
 *  @return unit part
 *  @see ex::content, ex::primpart, ex::unitcontprim */
ex ex::unit(const ex &x) const
{
	ex c = expand().lcoeff(x);
	if (is_exactly_a<numeric>(c))
		return c.info(info_flags::negative) ?_ex_1 : _ex1;
	else {
		ex y;
		if (get_first_symbol(c, y))
			return c.unit(y);
		else
			throw(std::invalid_argument("invalid expression in unit()"));
	}
}


/** Compute content part (= unit normal GCD of all coefficients) of a
 *  multivariate polynomial in Q[x]. The product of unit part, content part,
 *  and primitive part is the polynomial itself.
 *
 *  @param x  main variable
 *  @return content part
 *  @see ex::unit, ex::primpart, ex::unitcontprim */
ex ex::content(const ex &x) const
{
	if (is_exactly_a<numeric>(*this))
		return info(info_flags::negative) ? -*this : *this;

	ex e = expand();
	if (e.is_zero())
		return _ex0;

	// First, divide out the integer content (which we can calculate very efficiently).
	// If the leading coefficient of the quotient is an integer, we are done.
	ex c = e.integer_content();
	ex r = e / c;
	int deg = r.degree(x);
	ex lcoeff = r.coeff(x, deg);
	if (lcoeff.info(info_flags::integer))
		return c;

	// GCD of all coefficients
	int ldeg = r.ldegree(x);
	if (deg == ldeg)
		return lcoeff * c / lcoeff.unit(x);
	ex cont = _ex0;
	for (int i=ldeg; i<=deg; i++)
		cont = gcd(r.coeff(x, i), cont, NULL, NULL, false);
	return cont * c;
}


/** Compute primitive part of a multivariate polynomial in Q[x]. The result
 *  will be a unit-normal polynomial with a content part of 1. The product
 *  of unit part, content part, and primitive part is the polynomial itself.
 *
 *  @param x  main variable
 *  @return primitive part
 *  @see ex::unit, ex::content, ex::unitcontprim */
ex ex::primpart(const ex &x) const
{
	// We need to compute the unit and content anyway, so call unitcontprim()
	ex u, c, p;
	unitcontprim(x, u, c, p);
	return p;
}


/** Compute primitive part of a multivariate polynomial in Q[x] when the
 *  content part is already known. This function is faster in computing the
 *  primitive part than the previous function.
 *
 *  @param x  main variable
 *  @param c  previously computed content part
 *  @return primitive part */
ex ex::primpart(const ex &x, const ex &c) const
{
	if (is_zero() || c.is_zero())
		return _ex0;
	if (is_exactly_a<numeric>(*this))
		return _ex1;

	// Divide by unit and content to get primitive part
	ex u = unit(x);
	if (is_exactly_a<numeric>(c))
		return *this / (c * u);
	else
		return quo(*this, c * u, x, false);
}


/** Compute unit part, content part, and primitive part of a multivariate
 *  polynomial in Q[x]. The product of the three parts is the polynomial
 *  itself.
 *
 *  @param x  main variable
 *  @param u  unit part (returned)
 *  @param c  content part (returned)
 *  @param p  primitive part (returned)
 *  @see ex::unit, ex::content, ex::primpart */
void ex::unitcontprim(const ex &x, ex &u, ex &c, ex &p) const
{
	// Quick check for zero (avoid expanding)
	if (is_zero()) {
		u = _ex1;
		c = p = _ex0;
		return;
	}

	// Special case: input is a number
	if (is_exactly_a<numeric>(*this)) {
		if (info(info_flags::negative)) {
			u = _ex_1;
			c = abs(ex_to<numeric>(*this));
		} else {
			u = _ex1;
			c = *this;
		}
		p = _ex1;
		return;
	}

	// Expand input polynomial
	ex e = expand();
	if (e.is_zero()) {
		u = _ex1;
		c = p = _ex0;
		return;
	}

	// Compute unit and content
	u = unit(x);
	c = content(x);

	// Divide by unit and content to get primitive part
	if (c.is_zero()) {
		p = _ex0;
		return;
	}
	if (is_exactly_a<numeric>(c))
		p = *this / (c * u);
	else
		p = quo(e, c * u, x, false);
}


/*
 *  GCD of multivariate polynomials
 */

/** Compute GCD of multivariate polynomials using the subresultant PRS
 *  algorithm. This function is used internally by gcd().
 *
 *  @param a   first multivariate polynomial
 *  @param b   second multivariate polynomial
 *  @param var iterator to first element of vector of sym_desc structs
 *  @return the GCD as a new expression
 *  @see gcd */

static ex sr_gcd(const ex &a, const ex &b, sym_desc_vec::const_iterator var)
{
#if STATISTICS
	sr_gcd_called++;
#endif

	// The first symbol is our main variable
	const ex &x = var->sym;

	// Sort c and d so that c has higher degree
	ex c, d;
	int adeg = a.degree(x), bdeg = b.degree(x);
	int cdeg, ddeg;
	if (adeg >= bdeg) {
		c = a;
		d = b;
		cdeg = adeg;
		ddeg = bdeg;
	} else {
		c = b;
		d = a;
		cdeg = bdeg;
		ddeg = adeg;
	}

	// Remove content from c and d, to be attached to GCD later
	ex cont_c = c.content(x);
	ex cont_d = d.content(x);
	ex gamma = gcd(cont_c, cont_d, NULL, NULL, false);
	if (ddeg == 0)
		return gamma;
	c = c.primpart(x, cont_c);
	d = d.primpart(x, cont_d);

	// First element of subresultant sequence
	ex r = _ex0, ri = _ex1, psi = _ex1;
	int delta = cdeg - ddeg;

	for (;;) {

		// Calculate polynomial pseudo-remainder
		r = prem(c, d, x, false);
		if (r.is_zero())
			return gamma * d.primpart(x);

		c = d;
		cdeg = ddeg;
		if (!divide_in_z(r, ri * pow(psi, delta), d, var))
			throw(std::runtime_error("invalid expression in sr_gcd(), division failed"));
		ddeg = d.degree(x);
		if (ddeg == 0) {
			if (is_exactly_a<numeric>(r))
				return gamma;
			else
				return gamma * r.primpart(x);
		}

		// Next element of subresultant sequence
		ri = c.expand().lcoeff(x);
		if (delta == 1)
			psi = ri;
		else if (delta)
			divide_in_z(pow(ri, delta), pow(psi, delta-1), psi, var+1);
		delta = cdeg - ddeg;
	}
}


/** Return maximum (absolute value) coefficient of a polynomial.
 *  This function is used internally by heur_gcd().
 *
 *  @return maximum coefficient
 *  @see heur_gcd */
numeric ex::max_coefficient() const
{
	return bp->max_coefficient();
}

/** Implementation ex::max_coefficient().
 *  @see heur_gcd */
numeric basic::max_coefficient() const
{
	return *_num1_p;
}

numeric numeric::max_coefficient() const
{
	return abs(*this);
}

numeric add::max_coefficient() const
{
	epvector::const_iterator it = seq.begin();
	epvector::const_iterator itend = seq.end();
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	numeric cur_max = abs(ex_to<numeric>(overall_coeff));
	while (it != itend) {
		numeric a;
		GINAC_ASSERT(!is_exactly_a<numeric>(it->rest));
		a = abs(ex_to<numeric>(it->coeff));
		if (a > cur_max)
			cur_max = a;
		it++;
	}
	return cur_max;
}

numeric mul::max_coefficient() const
{
#ifdef DO_GINAC_ASSERT
	epvector::const_iterator it = seq.begin();
	epvector::const_iterator itend = seq.end();
	while (it != itend) {
		GINAC_ASSERT(!is_exactly_a<numeric>(recombine_pair_to_ex(*it)));
		it++;
	}
#endif // def DO_GINAC_ASSERT
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	return abs(ex_to<numeric>(overall_coeff));
}


/** Apply symmetric modular homomorphism to an expanded multivariate
 *  polynomial.  This function is usually used internally by heur_gcd().
 *
 *  @param xi  modulus
 *  @return mapped polynomial
 *  @see heur_gcd */
ex basic::smod(const numeric &xi) const
{
	return *this;
}

ex numeric::smod(const numeric &xi) const
{
	return GiNaC::smod(*this, xi);
}

ex add::smod(const numeric &xi) const
{
	epvector newseq;
	newseq.reserve(seq.size()+1);
	epvector::const_iterator it = seq.begin();
	epvector::const_iterator itend = seq.end();
	while (it != itend) {
		GINAC_ASSERT(!is_exactly_a<numeric>(it->rest));
		numeric coeff = GiNaC::smod(ex_to<numeric>(it->coeff), xi);
		if (!coeff.is_zero())
			newseq.push_back(expair(it->rest, coeff));
		it++;
	}
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	numeric coeff = GiNaC::smod(ex_to<numeric>(overall_coeff), xi);
	return (new add(newseq,coeff))->setflag(status_flags::dynallocated);
}

ex mul::smod(const numeric &xi) const
{
#ifdef DO_GINAC_ASSERT
	epvector::const_iterator it = seq.begin();
	epvector::const_iterator itend = seq.end();
	while (it != itend) {
		GINAC_ASSERT(!is_exactly_a<numeric>(recombine_pair_to_ex(*it)));
		it++;
	}
#endif // def DO_GINAC_ASSERT
	mul * mulcopyp = new mul(*this);
	GINAC_ASSERT(is_exactly_a<numeric>(overall_coeff));
	mulcopyp->overall_coeff = GiNaC::smod(ex_to<numeric>(overall_coeff),xi);
	mulcopyp->clearflag(status_flags::evaluated);
	mulcopyp->clearflag(status_flags::hash_calculated);
	return mulcopyp->setflag(status_flags::dynallocated);
}


/** xi-adic polynomial interpolation */
static ex interpolate(const ex &gamma, const numeric &xi, const ex &x, int degree_hint = 1)
{
	exvector g; g.reserve(degree_hint);
	ex e = gamma;
	numeric rxi = xi.inverse();
	for (int i=0; !e.is_zero(); i++) {
		ex gi = e.smod(xi);
		g.push_back(gi * power(x, i));
		e = (e - gi) * rxi;
	}
	return (new add(g))->setflag(status_flags::dynallocated);
}

/** Exception thrown by heur_gcd() to signal failure. */
class gcdheu_failed {};

/** Compute GCD of multivariate polynomials using the heuristic GCD algorithm.
 *  get_symbol_stats() must have been called previously with the input
 *  polynomials and an iterator to the first element of the sym_desc vector
 *  passed in. This function is used internally by gcd().
 *
 *  @param a  first multivariate polynomial (expanded)
 *  @param b  second multivariate polynomial (expanded)
 *  @param ca  cofactor of polynomial a (returned), NULL to suppress
 *             calculation of cofactor
 *  @param cb  cofactor of polynomial b (returned), NULL to suppress
 *             calculation of cofactor
 *  @param var iterator to first element of vector of sym_desc structs
 *  @return the GCD as a new expression
 *  @see gcd
 *  @exception gcdheu_failed() */
static ex heur_gcd(const ex &a, const ex &b, ex *ca, ex *cb, sym_desc_vec::const_iterator var)
{
#if STATISTICS
	heur_gcd_called++;
#endif

	// Algorithm only works for non-vanishing input polynomials
	if (a.is_zero() || b.is_zero())
		return (new fail())->setflag(status_flags::dynallocated);

	// GCD of two numeric values -> CLN
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b)) {
		numeric g = gcd(ex_to<numeric>(a), ex_to<numeric>(b));
		if (ca)
			*ca = ex_to<numeric>(a) / g;
		if (cb)
			*cb = ex_to<numeric>(b) / g;
		return g;
	}

	// The first symbol is our main variable
	const ex &x = var->sym;

	// Remove integer content
	numeric gc = gcd(a.integer_content(), b.integer_content());
	numeric rgc = gc.inverse();
	ex p = a * rgc;
	ex q = b * rgc;
	int maxdeg =  std::max(p.degree(x), q.degree(x));
	
	// Find evaluation point
	numeric mp = p.max_coefficient();
	numeric mq = q.max_coefficient();
	numeric xi;
	if (mp > mq)
		xi = mq * (*_num2_p) + (*_num2_p);
	else
		xi = mp * (*_num2_p) + (*_num2_p);

	// 6 tries maximum
	for (int t=0; t<6; t++) {
		if (xi.int_length() * maxdeg > 100000) {
			throw gcdheu_failed();
		}

		// Apply evaluation homomorphism and calculate GCD
		ex cp, cq;
		ex gamma = heur_gcd(p.subs(x == xi, subs_options::no_pattern), q.subs(x == xi, subs_options::no_pattern), &cp, &cq, var+1).expand();
		if (!is_exactly_a<fail>(gamma)) {

			// Reconstruct polynomial from GCD of mapped polynomials
			ex g = interpolate(gamma, xi, x, maxdeg);

			// Remove integer content
			g /= g.integer_content();

			// If the calculated polynomial divides both p and q, this is the GCD
			ex dummy;
			if (divide_in_z(p, g, ca ? *ca : dummy, var) && divide_in_z(q, g, cb ? *cb : dummy, var)) {
				g *= gc;
				return g;
			}
		}

		// Next evaluation point
		xi = iquo(xi * isqrt(isqrt(xi)) * numeric(73794), numeric(27011));
	}
	return (new fail())->setflag(status_flags::dynallocated);
}


/** Compute GCD (Greatest Common Divisor) of multivariate polynomials a(X)
 *  and b(X) in Z[X]. Optionally also compute the cofactors of a and b,
 *  defined by a = ca * gcd(a, b) and b = cb * gcd(a, b).
 *
 *  @param a  first multivariate polynomial
 *  @param b  second multivariate polynomial
 *  @param ca pointer to expression that will receive the cofactor of a, or NULL
 *  @param cb pointer to expression that will receive the cofactor of b, or NULL
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return the GCD as a new expression */
ex gcd(const ex &a, const ex &b, ex *ca, ex *cb, bool check_args)
{
#if STATISTICS
	gcd_called++;
#endif

	// GCD of numerics -> CLN
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b)) {
		numeric g = gcd(ex_to<numeric>(a), ex_to<numeric>(b));
		if (ca || cb) {
			if (g.is_zero()) {
				if (ca)
					*ca = _ex0;
				if (cb)
					*cb = _ex0;
			} else {
				if (ca)
					*ca = ex_to<numeric>(a) / g;
				if (cb)
					*cb = ex_to<numeric>(b) / g;
			}
		}
		return g;
	}

	// Check arguments
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial))) {
		throw(std::invalid_argument("gcd: arguments must be polynomials over the rationals"));
	}

	// Partially factored cases (to avoid expanding large expressions)
	if (is_exactly_a<mul>(a)) {
		if (is_exactly_a<mul>(b) && b.nops() > a.nops())
			goto factored_b;
factored_a:
		size_t num = a.nops();
		exvector g; g.reserve(num);
		exvector acc_ca; acc_ca.reserve(num);
		ex part_b = b;
		for (size_t i=0; i<num; i++) {
			ex part_ca, part_cb;
			g.push_back(gcd(a.op(i), part_b, &part_ca, &part_cb, check_args));
			acc_ca.push_back(part_ca);
			part_b = part_cb;
		}
		if (ca)
			*ca = (new mul(acc_ca))->setflag(status_flags::dynallocated);
		if (cb)
			*cb = part_b;
		return (new mul(g))->setflag(status_flags::dynallocated);
	} else if (is_exactly_a<mul>(b)) {
		if (is_exactly_a<mul>(a) && a.nops() > b.nops())
			goto factored_a;
factored_b:
		size_t num = b.nops();
		exvector g; g.reserve(num);
		exvector acc_cb; acc_cb.reserve(num);
		ex part_a = a;
		for (size_t i=0; i<num; i++) {
			ex part_ca, part_cb;
			g.push_back(gcd(part_a, b.op(i), &part_ca, &part_cb, check_args));
			acc_cb.push_back(part_cb);
			part_a = part_ca;
		}
		if (ca)
			*ca = part_a;
		if (cb)
			*cb = (new mul(acc_cb))->setflag(status_flags::dynallocated);
		return (new mul(g))->setflag(status_flags::dynallocated);
	}

#if FAST_COMPARE
	// Input polynomials of the form poly^n are sometimes also trivial
	if (is_exactly_a<power>(a)) {
		ex p = a.op(0);
		const ex& exp_a = a.op(1);
		if (is_exactly_a<power>(b)) {
			ex pb = b.op(0);
			const ex& exp_b = b.op(1);
			if (p.is_equal(pb)) {
				// a = p^n, b = p^m, gcd = p^min(n, m)
				if (exp_a < exp_b) {
					if (ca)
						*ca = _ex1;
					if (cb)
						*cb = power(p, exp_b - exp_a);
					return power(p, exp_a);
				} else {
					if (ca)
						*ca = power(p, exp_a - exp_b);
					if (cb)
						*cb = _ex1;
					return power(p, exp_b);
				}
			} else {
				ex p_co, pb_co;
				ex p_gcd = gcd(p, pb, &p_co, &pb_co, check_args);
				if (p_gcd.is_equal(_ex1)) {
					// a(x) = p(x)^n, b(x) = p_b(x)^m, gcd (p, p_b) = 1 ==>
					// gcd(a,b) = 1
					if (ca)
						*ca = a;
					if (cb)
						*cb = b;
					return _ex1;
					// XXX: do I need to check for p_gcd = -1?
				} else {
					// there are common factors:
					// a(x) = g(x)^n A(x)^n, b(x) = g(x)^m B(x)^m ==>
					// gcd(a, b) = g(x)^n gcd(A(x)^n, g(x)^(n-m) B(x)^m
					if (exp_a < exp_b) {
						return power(p_gcd, exp_a)*
							gcd(power(p_co, exp_a), power(p_gcd, exp_b-exp_a)*power(pb_co, exp_b), ca, cb, false);
					} else {
						return power(p_gcd, exp_b)*
							gcd(power(p_gcd, exp_a - exp_b)*power(p_co, exp_a), power(pb_co, exp_b), ca, cb, false);
					}
				} // p_gcd.is_equal(_ex1)
			} // p.is_equal(pb)

		} else {
			if (p.is_equal(b)) {
				// a = p^n, b = p, gcd = p
				if (ca)
					*ca = power(p, a.op(1) - 1);
				if (cb)
					*cb = _ex1;
				return p;
			} 

			ex p_co, bpart_co;
			ex p_gcd = gcd(p, b, &p_co, &bpart_co, false);

			if (p_gcd.is_equal(_ex1)) {
				// a(x) = p(x)^n, gcd(p, b) = 1 ==> gcd(a, b) = 1
				if (ca)
					*ca = a;
				if (cb)
					*cb = b;
				return _ex1;
			} else {
				// a(x) = g(x)^n A(x)^n, b(x) = g(x) B(x) ==> gcd(a, b) = g(x) gcd(g(x)^(n-1) A(x)^n, B(x))
				return p_gcd*gcd(power(p_gcd, exp_a-1)*power(p_co, exp_a), bpart_co, ca, cb, false);
			}
		} // is_exactly_a<power>(b)

	} else if (is_exactly_a<power>(b)) {
		ex p = b.op(0);
		if (p.is_equal(a)) {
			// a = p, b = p^n, gcd = p
			if (ca)
				*ca = _ex1;
			if (cb)
				*cb = power(p, b.op(1) - 1);
			return p;
		}

		ex p_co, apart_co;
		const ex& exp_b(b.op(1));
		ex p_gcd = gcd(a, p, &apart_co, &p_co, false);
		if (p_gcd.is_equal(_ex1)) {
			// b=p(x)^n, gcd(a, p) = 1 ==> gcd(a, b) == 1
			if (ca)
				*ca = a;
			if (cb)
				*cb = b;
			return _ex1;
		} else {
			// there are common factors:
			// a(x) = g(x) A(x), b(x) = g(x)^n B(x)^n ==> gcd = g(x) gcd(g(x)^(n-1) A(x)^n, B(x))

			return p_gcd*gcd(apart_co, power(p_gcd, exp_b-1)*power(p_co, exp_b), ca, cb, false);
		} // p_gcd.is_equal(_ex1)
	}
#endif

	// Some trivial cases
	ex aex = a.expand(), bex = b.expand();
	if (aex.is_zero()) {
		if (ca)
			*ca = _ex0;
		if (cb)
			*cb = _ex1;
		return b;
	}
	if (bex.is_zero()) {
		if (ca)
			*ca = _ex1;
		if (cb)
			*cb = _ex0;
		return a;
	}
	if (aex.is_equal(_ex1) || bex.is_equal(_ex1)) {
		if (ca)
			*ca = a;
		if (cb)
			*cb = b;
		return _ex1;
	}
#if FAST_COMPARE
	if (a.is_equal(b)) {
		if (ca)
			*ca = _ex1;
		if (cb)
			*cb = _ex1;
		return a;
	}
#endif

	if (is_a<symbol>(aex)) {
		if (! bex.subs(aex==_ex0, subs_options::no_pattern).is_zero()) {
			if (ca)
				*ca = a;
			if (cb)
				*cb = b;
			return _ex1;
		}
	}

	if (is_a<symbol>(bex)) {
		if (! aex.subs(bex==_ex0, subs_options::no_pattern).is_zero()) {
			if (ca)
				*ca = a;
			if (cb)
				*cb = b;
			return _ex1;
		}
	}

	if (is_exactly_a<numeric>(aex)) {
		numeric bcont = bex.integer_content();
		numeric g = gcd(ex_to<numeric>(aex), bcont);
		if (ca)
			*ca = ex_to<numeric>(aex)/g;
		if (cb)
	 		*cb = bex/g;
		return g;
	}

	if (is_exactly_a<numeric>(bex)) {
		numeric acont = aex.integer_content();
		numeric g = gcd(ex_to<numeric>(bex), acont);
		if (ca)
			*ca = aex/g;
		if (cb)
			*cb = ex_to<numeric>(bex)/g;
		return g;
	}

	// Gather symbol statistics
	sym_desc_vec sym_stats;
	get_symbol_stats(a, b, sym_stats);

	// The symbol with least degree which is contained in both polynomials
	// is our main variable
	sym_desc_vec::iterator vari = sym_stats.begin();
	while ((vari != sym_stats.end()) && 
	       (((vari->ldeg_b == 0) && (vari->deg_b == 0)) ||
	        ((vari->ldeg_a == 0) && (vari->deg_a == 0))))
		vari++;

	// No common symbols at all, just return 1:
	if (vari == sym_stats.end()) {
		// N.B: keep cofactors factored
		if (ca)
			*ca = a;
		if (cb)
			*cb = b;
		return _ex1;
	}
	// move symbols which contained only in one of the polynomials
	// to the end:
	rotate(sym_stats.begin(), vari, sym_stats.end());

	sym_desc_vec::const_iterator var = sym_stats.begin();
	const ex &x = var->sym;

	// Cancel trivial common factor
	int ldeg_a = var->ldeg_a;
	int ldeg_b = var->ldeg_b;
	int min_ldeg = std::min(ldeg_a,ldeg_b);
	if (min_ldeg > 0) {
		ex common = power(x, min_ldeg);
		return gcd((aex / common).expand(), (bex / common).expand(), ca, cb, false) * common;
	}

	// Try to eliminate variables
	if (var->deg_a == 0 && var->deg_b != 0 ) {
		ex bex_u, bex_c, bex_p;
		bex.unitcontprim(x, bex_u, bex_c, bex_p);
		ex g = gcd(aex, bex_c, ca, cb, false);
		if (cb)
			*cb *= bex_u * bex_p;
		return g;
	} else if (var->deg_b == 0 && var->deg_a != 0) {
		ex aex_u, aex_c, aex_p;
		aex.unitcontprim(x, aex_u, aex_c, aex_p);
		ex g = gcd(aex_c, bex, ca, cb, false);
		if (ca)
			*ca *= aex_u * aex_p;
		return g;
	}

	// Try heuristic algorithm first, fall back to PRS if that failed
	ex g;
	try {
		g = heur_gcd(aex, bex, ca, cb, var);
	} catch (gcdheu_failed) {
		g = fail();
	}
	if (is_exactly_a<fail>(g)) {
#if STATISTICS
		heur_gcd_failed++;
#endif
		g = sr_gcd(aex, bex, var);
		if (g.is_equal(_ex1)) {
			// Keep cofactors factored if possible
			if (ca)
				*ca = a;
			if (cb)
				*cb = b;
		} else {
			if (ca)
				divide(aex, g, *ca, false);
			if (cb)
				divide(bex, g, *cb, false);
		}
	} else {
		if (g.is_equal(_ex1)) {
			// Keep cofactors factored if possible
			if (ca)
				*ca = a;
			if (cb)
				*cb = b;
		}
	}

	return g;
}


/** Compute LCM (Least Common Multiple) of multivariate polynomials in Z[X].
 *
 *  @param a  first multivariate polynomial
 *  @param b  second multivariate polynomial
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return the LCM as a new expression */
ex lcm(const ex &a, const ex &b, bool check_args)
{
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b))
		return lcm(ex_to<numeric>(a), ex_to<numeric>(b));
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("lcm: arguments must be polynomials over the rationals"));
	
	ex ca, cb;
	ex g = gcd(a, b, &ca, &cb, false);
	return ca * cb * g;
}


/*
 *  Square-free factorization
 */

/** Compute square-free factorization of multivariate polynomial a(x) using
 *  Yun's algorithm.  Used internally by sqrfree().
 *
 *  @param a  multivariate polynomial over Z[X], treated here as univariate
 *            polynomial in x.
 *  @param x  variable to factor in
 *  @return   vector of factors sorted in ascending degree */
static exvector sqrfree_yun(const ex &a, const symbol &x)
{
	exvector res;
	ex w = a;
	ex z = w.diff(x);
	ex g = gcd(w, z);
	if (g.is_equal(_ex1)) {
		res.push_back(a);
		return res;
	}
	ex y;
	do {
		w = quo(w, g, x);
		y = quo(z, g, x);
		z = y - w.diff(x);
		g = gcd(w, z);
		res.push_back(g);
	} while (!z.is_zero());
	return res;
}


/** Compute a square-free factorization of a multivariate polynomial in Q[X].
 *
 *  @param a  multivariate polynomial over Q[X]
 *  @param l  lst of variables to factor in, may be left empty for autodetection
 *  @return   a square-free factorization of \p a.
 *
 * \note
 * A polynomial \f$p(X) \in C[X]\f$ is said <EM>square-free</EM>
 * if, whenever any two polynomials \f$q(X)\f$ and \f$r(X)\f$
 * are such that
 * \f[
 *     p(X) = q(X)^2 r(X),
 * \f]
 * we have \f$q(X) \in C\f$.
 * This means that \f$p(X)\f$ has no repeated factors, apart
 * eventually from constants.
 * Given a polynomial \f$p(X) \in C[X]\f$, we say that the
 * decomposition
 * \f[
 *   p(X) = b \cdot p_1(X)^{a_1} \cdot p_2(X)^{a_2} \cdots p_r(X)^{a_r}
 * \f]
 * is a <EM>square-free factorization</EM> of \f$p(X)\f$ if the
 * following conditions hold:
 * -#  \f$b \in C\f$ and \f$b \neq 0\f$;
 * -#  \f$a_i\f$ is a positive integer for \f$i = 1, \ldots, r\f$;
 * -#  the degree of the polynomial \f$p_i\f$ is strictly positive
 *     for \f$i = 1, \ldots, r\f$;
 * -#  the polynomial \f$\Pi_{i=1}^r p_i(X)\f$ is square-free.
 *
 * Square-free factorizations need not be unique.  For example, if
 * \f$a_i\f$ is even, we could change the polynomial \f$p_i(X)\f$
 * into \f$-p_i(X)\f$.
 * Observe also that the factors \f$p_i(X)\f$ need not be irreducible
 * polynomials.
 */
ex sqrfree(const ex &a, const lst &l)
{
	if (is_exactly_a<numeric>(a) ||     // algorithm does not trap a==0
	    is_a<symbol>(a))        // shortcut
		return a;

	// If no lst of variables to factorize in was specified we have to
	// invent one now.  Maybe one can optimize here by reversing the order
	// or so, I don't know.
	lst args;
	if (l.nops()==0) {
		sym_desc_vec sdv;
		get_symbol_stats(a, _ex0, sdv);
		sym_desc_vec::const_iterator it = sdv.begin(), itend = sdv.end();
		while (it != itend) {
			args.append(it->sym);
			++it;
		}
	} else {
		args = l;
	}

	// Find the symbol to factor in at this stage
	if (!is_a<symbol>(args.op(0)))
		throw (std::runtime_error("sqrfree(): invalid factorization variable"));
	const symbol &x = ex_to<symbol>(args.op(0));

	// convert the argument from something in Q[X] to something in Z[X]
	const numeric lcm = lcm_of_coefficients_denominators(a);
	const ex tmp = multiply_lcm(a,lcm);

	// find the factors
	exvector factors = sqrfree_yun(tmp, x);

	// construct the next list of symbols with the first element popped
	lst newargs = args;
	newargs.remove_first();

	// recurse down the factors in remaining variables
	if (newargs.nops()>0) {
		exvector::iterator i = factors.begin();
		while (i != factors.end()) {
			*i = sqrfree(*i, newargs);
			++i;
		}
	}

	// Done with recursion, now construct the final result
	ex result = _ex1;
	exvector::const_iterator it = factors.begin(), itend = factors.end();
	for (int p = 1; it!=itend; ++it, ++p)
		result *= power(*it, p);

	// Yun's algorithm does not account for constant factors.  (For univariate
	// polynomials it works only in the monic case.)  We can correct this by
	// inserting what has been lost back into the result.  For completeness
	// we'll also have to recurse down that factor in the remaining variables.
	if (newargs.nops()>0)
		result *= sqrfree(quo(tmp, result, x), newargs);
	else
		result *= quo(tmp, result, x);

	// Put in the reational overall factor again and return
	return result *	lcm.inverse();
}


/** Compute square-free partial fraction decomposition of rational function
 *  a(x).
 *
 *  @param a rational function over Z[x], treated as univariate polynomial
 *           in x
 *  @param x variable to factor in
 *  @return decomposed rational function */
ex sqrfree_parfrac(const ex & a, const symbol & x)
{
	// Find numerator and denominator
	ex nd = numer_denom(a);
	ex numer = nd.op(0), denom = nd.op(1);
//clog << "numer = " << numer << ", denom = " << denom << endl;

	// Convert N(x)/D(x) -> Q(x) + R(x)/D(x), so degree(R) < degree(D)
	ex red_poly = quo(numer, denom, x), red_numer = rem(numer, denom, x).expand();
//clog << "red_poly = " << red_poly << ", red_numer = " << red_numer << endl;

	// Factorize denominator and compute cofactors
	exvector yun = sqrfree_yun(denom, x);
//clog << "yun factors: " << exprseq(yun) << endl;
	size_t num_yun = yun.size();
	exvector factor; factor.reserve(num_yun);
	exvector cofac; cofac.reserve(num_yun);
	for (size_t i=0; i<num_yun; i++) {
		if (!yun[i].is_equal(_ex1)) {
			for (size_t j=0; j<=i; j++) {
				factor.push_back(pow(yun[i], j+1));
				ex prod = _ex1;
				for (size_t k=0; k<num_yun; k++) {
					if (k == i)
						prod *= pow(yun[k], i-j);
					else
						prod *= pow(yun[k], k+1);
				}
				cofac.push_back(prod.expand());
			}
		}
	}
	size_t num_factors = factor.size();
//clog << "factors  : " << exprseq(factor) << endl;
//clog << "cofactors: " << exprseq(cofac) << endl;

	// Construct coefficient matrix for decomposition
	int max_denom_deg = denom.degree(x);
	matrix sys(max_denom_deg + 1, num_factors);
	matrix rhs(max_denom_deg + 1, 1);
	for (int i=0; i<=max_denom_deg; i++) {
		for (size_t j=0; j<num_factors; j++)
			sys(i, j) = cofac[j].coeff(x, i);
		rhs(i, 0) = red_numer.coeff(x, i);
	}
//clog << "coeffs: " << sys << endl;
//clog << "rhs   : " << rhs << endl;

	// Solve resulting linear system
	matrix vars(num_factors, 1);
	for (size_t i=0; i<num_factors; i++)
		vars(i, 0) = symbol();
	matrix sol = sys.solve(vars, rhs);

	// Sum up decomposed fractions
	ex sum = 0;
	for (size_t i=0; i<num_factors; i++)
		sum += sol(i, 0) / factor[i];

	return red_poly + sum;
}


/*
 *  Normal form of rational functions
 */

/*
 *  Note: The internal normal() functions (= basic::normal() and overloaded
 *  functions) all return lists of the form {numerator, denominator}. This
 *  is to get around mul::eval()'s automatic expansion of numeric coefficients.
 *  E.g. (a+b)/3 is automatically converted to a/3+b/3 but we want to keep
 *  the information that (a+b) is the numerator and 3 is the denominator.
 */


/** Create a symbol for replacing the expression "e" (or return a previously
 *  assigned symbol). The symbol and expression are appended to repl, for
 *  a later application of subs().
 *  @see ex::normal */
static ex replace_with_symbol(const ex & e, exmap & repl, exmap & rev_lookup)
{
	// Expression already replaced? Then return the assigned symbol
	exmap::const_iterator it = rev_lookup.find(e);
	if (it != rev_lookup.end())
		return it->second;
	
	// Otherwise create new symbol and add to list, taking care that the
	// replacement expression doesn't itself contain symbols from repl,
	// because subs() is not recursive
	ex es = (new symbol)->setflag(status_flags::dynallocated);
	ex e_replaced = e.subs(repl, subs_options::no_pattern);
	repl.insert(std::make_pair(es, e_replaced));
	rev_lookup.insert(std::make_pair(e_replaced, es));
	return es;
}

/** Create a symbol for replacing the expression "e" (or return a previously
 *  assigned symbol). The symbol and expression are appended to repl, and the
 *  symbol is returned.
 *  @see basic::to_rational
 *  @see basic::to_polynomial */
static ex replace_with_symbol(const ex & e, exmap & repl)
{
	// Expression already replaced? Then return the assigned symbol
	for (exmap::const_iterator it = repl.begin(); it != repl.end(); ++it)
		if (it->second.is_equal(e))
			return it->first;
	
	// Otherwise create new symbol and add to list, taking care that the
	// replacement expression doesn't itself contain symbols from repl,
	// because subs() is not recursive
	ex es = (new symbol)->setflag(status_flags::dynallocated);
	ex e_replaced = e.subs(repl, subs_options::no_pattern);
	repl.insert(std::make_pair(es, e_replaced));
	return es;
}


/** Function object to be applied by basic::normal(). */
struct normal_map_function : public map_function {
	int level;
	normal_map_function(int l) : level(l) {}
	ex operator()(const ex & e) { return normal(e, level); }
};

/** Default implementation of ex::normal(). It normalizes the children and
 *  replaces the object with a temporary symbol.
 *  @see ex::normal */
ex basic::normal(exmap & repl, exmap & rev_lookup, int level) const
{
	if (nops() == 0)
		return (new lst(replace_with_symbol(*this, repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
	else {
		if (level == 1)
			return (new lst(replace_with_symbol(*this, repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
		else if (level == -max_recursion_level)
			throw(std::runtime_error("max recursion level reached"));
		else {
			normal_map_function map_normal(level - 1);
			return (new lst(replace_with_symbol(map(map_normal), repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
		}
	}
}


/** Implementation of ex::normal() for symbols. This returns the unmodified symbol.
 *  @see ex::normal */
ex symbol::normal(exmap & repl, exmap & rev_lookup, int level) const
{
	return (new lst(*this, _ex1))->setflag(status_flags::dynallocated);
}


/** Implementation of ex::normal() for a numeric. It splits complex numbers
 *  into re+I*im and replaces I and non-rational real numbers with a temporary
 *  symbol.
 *  @see ex::normal */
ex numeric::normal(exmap & repl, exmap & rev_lookup, int level) const
{
	numeric num = numer();
	ex numex = num;

	if (num.is_real()) {
		if (!num.is_integer())
			numex = replace_with_symbol(numex, repl, rev_lookup);
	} else { // complex
		numeric re = num.real(), im = num.imag();
		ex re_ex = re.is_rational() ? re : replace_with_symbol(re, repl, rev_lookup);
		ex im_ex = im.is_rational() ? im : replace_with_symbol(im, repl, rev_lookup);
		numex = re_ex + im_ex * replace_with_symbol(I, repl, rev_lookup);
	}

	// Denominator is always a real integer (see numeric::denom())
	return (new lst(numex, denom()))->setflag(status_flags::dynallocated);
}


/** Fraction cancellation.
 *  @param n  numerator
 *  @param d  denominator
 *  @return cancelled fraction {n, d} as a list */
static ex frac_cancel(const ex &n, const ex &d)
{
	ex num = n;
	ex den = d;
	numeric pre_factor = *_num1_p;

//std::clog << "frac_cancel num = " << num << ", den = " << den << std::endl;

	// Handle trivial case where denominator is 1
	if (den.is_equal(_ex1))
		return (new lst(num, den))->setflag(status_flags::dynallocated);

	// Handle special cases where numerator or denominator is 0
	if (num.is_zero())
		return (new lst(num, _ex1))->setflag(status_flags::dynallocated);
	if (den.expand().is_zero())
		throw(std::overflow_error("frac_cancel: division by zero in frac_cancel"));

	// Bring numerator and denominator to Z[X] by multiplying with
	// LCM of all coefficients' denominators
	numeric num_lcm = lcm_of_coefficients_denominators(num);
	numeric den_lcm = lcm_of_coefficients_denominators(den);
	num = multiply_lcm(num, num_lcm);
	den = multiply_lcm(den, den_lcm);
	pre_factor = den_lcm / num_lcm;

	// Cancel GCD from numerator and denominator
	ex cnum, cden;
	if (gcd(num, den, &cnum, &cden, false) != _ex1) {
		num = cnum;
		den = cden;
	}

	// Make denominator unit normal (i.e. coefficient of first symbol
	// as defined by get_first_symbol() is made positive)
	if (is_exactly_a<numeric>(den)) {
		if (ex_to<numeric>(den).is_negative()) {
			num *= _ex_1;
			den *= _ex_1;
		}
	} else {
		ex x;
		if (get_first_symbol(den, x)) {
			GINAC_ASSERT(is_exactly_a<numeric>(den.unit(x)));
			if (ex_to<numeric>(den.unit(x)).is_negative()) {
				num *= _ex_1;
				den *= _ex_1;
			}
		}
	}

	// Return result as list
//std::clog << " returns num = " << num << ", den = " << den << ", pre_factor = " << pre_factor << std::endl;
	return (new lst(num * pre_factor.numer(), den * pre_factor.denom()))->setflag(status_flags::dynallocated);
}


/** Implementation of ex::normal() for a sum. It expands terms and performs
 *  fractional addition.
 *  @see ex::normal */
ex add::normal(exmap & repl, exmap & rev_lookup, int level) const
{
	if (level == 1)
		return (new lst(replace_with_symbol(*this, repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
	else if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));

	// Normalize children and split each one into numerator and denominator
	exvector nums, dens;
	nums.reserve(seq.size()+1);
	dens.reserve(seq.size()+1);
	epvector::const_iterator it = seq.begin(), itend = seq.end();
	while (it != itend) {
		ex n = ex_to<basic>(recombine_pair_to_ex(*it)).normal(repl, rev_lookup, level-1);
		nums.push_back(n.op(0));
		dens.push_back(n.op(1));
		it++;
	}
	ex n = ex_to<numeric>(overall_coeff).normal(repl, rev_lookup, level-1);
	nums.push_back(n.op(0));
	dens.push_back(n.op(1));
	GINAC_ASSERT(nums.size() == dens.size());

	// Now, nums is a vector of all numerators and dens is a vector of
	// all denominators
//std::clog << "add::normal uses " << nums.size() << " summands:\n";

	// Add fractions sequentially
	exvector::const_iterator num_it = nums.begin(), num_itend = nums.end();
	exvector::const_iterator den_it = dens.begin(), den_itend = dens.end();
//std::clog << " num = " << *num_it << ", den = " << *den_it << std::endl;
	ex num = *num_it++, den = *den_it++;
	while (num_it != num_itend) {
//std::clog << " num = " << *num_it << ", den = " << *den_it << std::endl;
		ex next_num = *num_it++, next_den = *den_it++;

		// Trivially add sequences of fractions with identical denominators
		while ((den_it != den_itend) && next_den.is_equal(*den_it)) {
			next_num += *num_it;
			num_it++; den_it++;
		}

		// Additiion of two fractions, taking advantage of the fact that
		// the heuristic GCD algorithm computes the cofactors at no extra cost
		ex co_den1, co_den2;
		ex g = gcd(den, next_den, &co_den1, &co_den2, false);
		num = ((num * co_den2) + (next_num * co_den1)).expand();
		den *= co_den2;		// this is the lcm(den, next_den)
	}
//std::clog << " common denominator = " << den << std::endl;

	// Cancel common factors from num/den
	return frac_cancel(num, den);
}


/** Implementation of ex::normal() for a product. It cancels common factors
 *  from fractions.
 *  @see ex::normal() */
ex mul::normal(exmap & repl, exmap & rev_lookup, int level) const
{
	if (level == 1)
		return (new lst(replace_with_symbol(*this, repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
	else if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));

	// Normalize children, separate into numerator and denominator
	exvector num; num.reserve(seq.size());
	exvector den; den.reserve(seq.size());
	ex n;
	epvector::const_iterator it = seq.begin(), itend = seq.end();
	while (it != itend) {
		n = ex_to<basic>(recombine_pair_to_ex(*it)).normal(repl, rev_lookup, level-1);
		num.push_back(n.op(0));
		den.push_back(n.op(1));
		it++;
	}
	n = ex_to<numeric>(overall_coeff).normal(repl, rev_lookup, level-1);
	num.push_back(n.op(0));
	den.push_back(n.op(1));

	// Perform fraction cancellation
	return frac_cancel((new mul(num))->setflag(status_flags::dynallocated),
	                   (new mul(den))->setflag(status_flags::dynallocated));
}


/** Implementation of ex::normal([B) for powers. It normalizes the basis,
 *  distributes integer exponents to numerator and denominator, and replaces
 *  non-integer powers by temporary symbols.
 *  @see ex::normal */
ex power::normal(exmap & repl, exmap & rev_lookup, int level) const
{
	if (level == 1)
		return (new lst(replace_with_symbol(*this, repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
	else if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));

	// Normalize basis and exponent (exponent gets reassembled)
	ex n_basis = ex_to<basic>(basis).normal(repl, rev_lookup, level-1);
	ex n_exponent = ex_to<basic>(exponent).normal(repl, rev_lookup, level-1);
	n_exponent = n_exponent.op(0) / n_exponent.op(1);

	if (n_exponent.info(info_flags::integer)) {

		if (n_exponent.info(info_flags::positive)) {

			// (a/b)^n -> {a^n, b^n}
			return (new lst(power(n_basis.op(0), n_exponent), power(n_basis.op(1), n_exponent)))->setflag(status_flags::dynallocated);

		} else if (n_exponent.info(info_flags::negative)) {

			// (a/b)^-n -> {b^n, a^n}
			return (new lst(power(n_basis.op(1), -n_exponent), power(n_basis.op(0), -n_exponent)))->setflag(status_flags::dynallocated);
		}

	} else {

		if (n_exponent.info(info_flags::positive)) {

			// (a/b)^x -> {sym((a/b)^x), 1}
			return (new lst(replace_with_symbol(power(n_basis.op(0) / n_basis.op(1), n_exponent), repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);

		} else if (n_exponent.info(info_flags::negative)) {

			if (n_basis.op(1).is_equal(_ex1)) {

				// a^-x -> {1, sym(a^x)}
				return (new lst(_ex1, replace_with_symbol(power(n_basis.op(0), -n_exponent), repl, rev_lookup)))->setflag(status_flags::dynallocated);

			} else {

				// (a/b)^-x -> {sym((b/a)^x), 1}
				return (new lst(replace_with_symbol(power(n_basis.op(1) / n_basis.op(0), -n_exponent), repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
			}
		}
	}

	// (a/b)^x -> {sym((a/b)^x, 1}
	return (new lst(replace_with_symbol(power(n_basis.op(0) / n_basis.op(1), n_exponent), repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
}


/** Implementation of ex::normal() for pseries. It normalizes each coefficient
 *  and replaces the series by a temporary symbol.
 *  @see ex::normal */
ex pseries::normal(exmap & repl, exmap & rev_lookup, int level) const
{
	epvector newseq;
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		ex restexp = i->rest.normal();
		if (!restexp.is_zero())
			newseq.push_back(expair(restexp, i->coeff));
		++i;
	}
	ex n = pseries(relational(var,point), newseq);
	return (new lst(replace_with_symbol(n, repl, rev_lookup), _ex1))->setflag(status_flags::dynallocated);
}


/** Normalization of rational functions.
 *  This function converts an expression to its normal form
 *  "numerator/denominator", where numerator and denominator are (relatively
 *  prime) polynomials. Any subexpressions which are not rational functions
 *  (like non-rational numbers, non-integer powers or functions like sin(),
 *  cos() etc.) are replaced by temporary symbols which are re-substituted by
 *  the (normalized) subexpressions before normal() returns (this way, any
 *  expression can be treated as a rational function). normal() is applied
 *  recursively to arguments of functions etc.
 *
 *  @param level maximum depth of recursion
 *  @return normalized expression */
ex ex::normal(int level) const
{
	exmap repl, rev_lookup;

	ex e = bp->normal(repl, rev_lookup, level);
	GINAC_ASSERT(is_a<lst>(e));

	// Re-insert replaced symbols
	if (!repl.empty())
		e = e.subs(repl, subs_options::no_pattern);

	// Convert {numerator, denominator} form back to fraction
	return e.op(0) / e.op(1);
}

/** Get numerator of an expression. If the expression is not of the normal
 *  form "numerator/denominator", it is first converted to this form and
 *  then the numerator is returned.
 *
 *  @see ex::normal
 *  @return numerator */
ex ex::numer() const
{
	exmap repl, rev_lookup;

	ex e = bp->normal(repl, rev_lookup, 0);
	GINAC_ASSERT(is_a<lst>(e));

	// Re-insert replaced symbols
	if (repl.empty())
		return e.op(0);
	else
		return e.op(0).subs(repl, subs_options::no_pattern);
}

/** Get denominator of an expression. If the expression is not of the normal
 *  form "numerator/denominator", it is first converted to this form and
 *  then the denominator is returned.
 *
 *  @see ex::normal
 *  @return denominator */
ex ex::denom() const
{
	exmap repl, rev_lookup;

	ex e = bp->normal(repl, rev_lookup, 0);
	GINAC_ASSERT(is_a<lst>(e));

	// Re-insert replaced symbols
	if (repl.empty())
		return e.op(1);
	else
		return e.op(1).subs(repl, subs_options::no_pattern);
}

/** Get numerator and denominator of an expression. If the expresison is not
 *  of the normal form "numerator/denominator", it is first converted to this
 *  form and then a list [numerator, denominator] is returned.
 *
 *  @see ex::normal
 *  @return a list [numerator, denominator] */
ex ex::numer_denom() const
{
	exmap repl, rev_lookup;

	ex e = bp->normal(repl, rev_lookup, 0);
	GINAC_ASSERT(is_a<lst>(e));

	// Re-insert replaced symbols
	if (repl.empty())
		return e;
	else
		return e.subs(repl, subs_options::no_pattern);
}


/** Rationalization of non-rational functions.
 *  This function converts a general expression to a rational function
 *  by replacing all non-rational subexpressions (like non-rational numbers,
 *  non-integer powers or functions like sin(), cos() etc.) to temporary
 *  symbols. This makes it possible to use functions like gcd() and divide()
 *  on non-rational functions by applying to_rational() on the arguments,
 *  calling the desired function and re-substituting the temporary symbols
 *  in the result. To make the last step possible, all temporary symbols and
 *  their associated expressions are collected in the map specified by the
 *  repl parameter, ready to be passed as an argument to ex::subs().
 *
 *  @param repl collects all temporary symbols and their replacements
 *  @return rationalized expression */
ex ex::to_rational(exmap & repl) const
{
	return bp->to_rational(repl);
}

// GiNaC 1.1 compatibility function
ex ex::to_rational(lst & repl_lst) const
{
	// Convert lst to exmap
	exmap m;
	for (lst::const_iterator it = repl_lst.begin(); it != repl_lst.end(); ++it)
		m.insert(std::make_pair(it->op(0), it->op(1)));

	ex ret = bp->to_rational(m);

	// Convert exmap back to lst
	repl_lst.remove_all();
	for (exmap::const_iterator it = m.begin(); it != m.end(); ++it)
		repl_lst.append(it->first == it->second);

	return ret;
}

ex ex::to_polynomial(exmap & repl) const
{
	return bp->to_polynomial(repl);
}

// GiNaC 1.1 compatibility function
ex ex::to_polynomial(lst & repl_lst) const
{
	// Convert lst to exmap
	exmap m;
	for (lst::const_iterator it = repl_lst.begin(); it != repl_lst.end(); ++it)
		m.insert(std::make_pair(it->op(0), it->op(1)));

	ex ret = bp->to_polynomial(m);

	// Convert exmap back to lst
	repl_lst.remove_all();
	for (exmap::const_iterator it = m.begin(); it != m.end(); ++it)
		repl_lst.append(it->first == it->second);

	return ret;
}

/** Default implementation of ex::to_rational(). This replaces the object with
 *  a temporary symbol. */
ex basic::to_rational(exmap & repl) const
{
	return replace_with_symbol(*this, repl);
}

ex basic::to_polynomial(exmap & repl) const
{
	return replace_with_symbol(*this, repl);
}


/** Implementation of ex::to_rational() for symbols. This returns the
 *  unmodified symbol. */
ex symbol::to_rational(exmap & repl) const
{
	return *this;
}

/** Implementation of ex::to_polynomial() for symbols. This returns the
 *  unmodified symbol. */
ex symbol::to_polynomial(exmap & repl) const
{
	return *this;
}


/** Implementation of ex::to_rational() for a numeric. It splits complex
 *  numbers into re+I*im and replaces I and non-rational real numbers with a
 *  temporary symbol. */
ex numeric::to_rational(exmap & repl) const
{
	if (is_real()) {
		if (!is_rational())
			return replace_with_symbol(*this, repl);
	} else { // complex
		numeric re = real();
		numeric im = imag();
		ex re_ex = re.is_rational() ? re : replace_with_symbol(re, repl);
		ex im_ex = im.is_rational() ? im : replace_with_symbol(im, repl);
		return re_ex + im_ex * replace_with_symbol(I, repl);
	}
	return *this;
}

/** Implementation of ex::to_polynomial() for a numeric. It splits complex
 *  numbers into re+I*im and replaces I and non-integer real numbers with a
 *  temporary symbol. */
ex numeric::to_polynomial(exmap & repl) const
{
	if (is_real()) {
		if (!is_integer())
			return replace_with_symbol(*this, repl);
	} else { // complex
		numeric re = real();
		numeric im = imag();
		ex re_ex = re.is_integer() ? re : replace_with_symbol(re, repl);
		ex im_ex = im.is_integer() ? im : replace_with_symbol(im, repl);
		return re_ex + im_ex * replace_with_symbol(I, repl);
	}
	return *this;
}


/** Implementation of ex::to_rational() for powers. It replaces non-integer
 *  powers by temporary symbols. */
ex power::to_rational(exmap & repl) const
{
	if (exponent.info(info_flags::integer))
		return power(basis.to_rational(repl), exponent);
	else
		return replace_with_symbol(*this, repl);
}

/** Implementation of ex::to_polynomial() for powers. It replaces non-posint
 *  powers by temporary symbols. */
ex power::to_polynomial(exmap & repl) const
{
	if (exponent.info(info_flags::posint))
		return power(basis.to_rational(repl), exponent);
	else if (exponent.info(info_flags::negint))
	{
		ex basis_pref = collect_common_factors(basis);
		if (is_exactly_a<mul>(basis_pref) || is_exactly_a<power>(basis_pref)) {
			// (A*B)^n will be automagically transformed to A^n*B^n
			ex t = power(basis_pref, exponent);
			return t.to_polynomial(repl);
		}
		else
			return power(replace_with_symbol(power(basis, _ex_1), repl), -exponent);
	} 
	else
		return replace_with_symbol(*this, repl);
}


/** Implementation of ex::to_rational() for expairseqs. */
ex expairseq::to_rational(exmap & repl) const
{
	epvector s;
	s.reserve(seq.size());
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		s.push_back(split_ex_to_pair(recombine_pair_to_ex(*i).to_rational(repl)));
		++i;
	}
	ex oc = overall_coeff.to_rational(repl);
	if (oc.info(info_flags::numeric))
		return thisexpairseq(s, overall_coeff);
	else
		s.push_back(combine_ex_with_coeff_to_pair(oc, _ex1));
	return thisexpairseq(s, default_overall_coeff());
}

/** Implementation of ex::to_polynomial() for expairseqs. */
ex expairseq::to_polynomial(exmap & repl) const
{
	epvector s;
	s.reserve(seq.size());
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		s.push_back(split_ex_to_pair(recombine_pair_to_ex(*i).to_polynomial(repl)));
		++i;
	}
	ex oc = overall_coeff.to_polynomial(repl);
	if (oc.info(info_flags::numeric))
		return thisexpairseq(s, overall_coeff);
	else
		s.push_back(combine_ex_with_coeff_to_pair(oc, _ex1));
	return thisexpairseq(s, default_overall_coeff());
}


/** Remove the common factor in the terms of a sum 'e' by calculating the GCD,
 *  and multiply it into the expression 'factor' (which needs to be initialized
 *  to 1, unless you're accumulating factors). */
static ex find_common_factor(const ex & e, ex & factor, exmap & repl)
{
	if (is_exactly_a<add>(e)) {

		size_t num = e.nops();
		exvector terms; terms.reserve(num);
		ex gc;

		// Find the common GCD
		for (size_t i=0; i<num; i++) {
			ex x = e.op(i).to_polynomial(repl);

			if (is_exactly_a<add>(x) || is_exactly_a<mul>(x) || is_a<power>(x)) {
				ex f = 1;
				x = find_common_factor(x, f, repl);
				x *= f;
			}

			if (i == 0)
				gc = x;
			else
				gc = gcd(gc, x);

			terms.push_back(x);
		}

		if (gc.is_equal(_ex1))
			return e;

		// The GCD is the factor we pull out
		factor *= gc;

		// Now divide all terms by the GCD
		for (size_t i=0; i<num; i++) {
			ex x;

			// Try to avoid divide() because it expands the polynomial
			ex &t = terms[i];
			if (is_exactly_a<mul>(t)) {
				for (size_t j=0; j<t.nops(); j++) {
					if (t.op(j).is_equal(gc)) {
						exvector v; v.reserve(t.nops());
						for (size_t k=0; k<t.nops(); k++) {
							if (k == j)
								v.push_back(_ex1);
							else
								v.push_back(t.op(k));
						}
						t = (new mul(v))->setflag(status_flags::dynallocated);
						goto term_done;
					}
				}
			}

			divide(t, gc, x);
			t = x;
term_done:	;
		}
		return (new add(terms))->setflag(status_flags::dynallocated);

	} else if (is_exactly_a<mul>(e)) {

		size_t num = e.nops();
		exvector v; v.reserve(num);

		for (size_t i=0; i<num; i++)
			v.push_back(find_common_factor(e.op(i), factor, repl));

		return (new mul(v))->setflag(status_flags::dynallocated);

	} else if (is_exactly_a<power>(e)) {
		const ex e_exp(e.op(1));
		if (e_exp.info(info_flags::integer)) {
			ex eb = e.op(0).to_polynomial(repl);
			ex factor_local(_ex1);
			ex pre_res = find_common_factor(eb, factor_local, repl);
			factor *= power(factor_local, e_exp);
			return power(pre_res, e_exp);
			
		} else
			return e.to_polynomial(repl);

	} else
		return e;
}


/** Collect common factors in sums. This converts expressions like
 *  'a*(b*x+b*y)' to 'a*b*(x+y)'. */
ex collect_common_factors(const ex & e)
{
	if (is_exactly_a<add>(e) || is_exactly_a<mul>(e) || is_exactly_a<power>(e)) {

		exmap repl;
		ex factor = 1;
		ex r = find_common_factor(e, factor, repl);
		return factor.subs(repl, subs_options::no_pattern) * r.subs(repl, subs_options::no_pattern);

	} else
		return e;
}


/** Resultant of two expressions e1,e2 with respect to symbol s.
 *  Method: Compute determinant of Sylvester matrix of e1,e2,s.  */
ex resultant(const ex & e1, const ex & e2, const ex & s)
{
	const ex ee1 = e1.expand();
	const ex ee2 = e2.expand();
	if (!ee1.info(info_flags::polynomial) ||
	    !ee2.info(info_flags::polynomial))
		throw(std::runtime_error("resultant(): arguments must be polynomials"));

	const int h1 = ee1.degree(s);
	const int l1 = ee1.ldegree(s);
	const int h2 = ee2.degree(s);
	const int l2 = ee2.ldegree(s);

	const int msize = h1 + h2;
	matrix m(msize, msize);

	for (int l = h1; l >= l1; --l) {
		const ex e = ee1.coeff(s, l);
		for (int k = 0; k < h2; ++k)
			m(k, k+h1-l) = e;
	}
	for (int l = h2; l >= l2; --l) {
		const ex e = ee2.coeff(s, l);
		for (int k = 0; k < h1; ++k)
			m(k+h2, k+h2-l) = e;
	}

	return m.determinant();
}


} // namespace GiNaC
