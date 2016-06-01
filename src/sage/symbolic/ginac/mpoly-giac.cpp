/** @file mpoly-giac.cpp
 * 
 *  Copyright (C) 2016  Ralf Stephan
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

#ifdef HAVE_CONFIG_H
#include "pynac-config.h"
#endif

#ifdef PYNAC_HAVE_LIBGIAC

#include <string>
#include <iostream>
#include <sstream>
#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE

#include "upoly.h"
#include "basic.h"
#include "ex.h"
#include "add.h"
#include "constant.h"
#include "expairseq.h"
#include "mul.h"
#include "numeric.h"
#include "power.h"
#include "operators.h"
#include "pseries.h"
#include "symbol.h"
#include "function.h"
#include "utils.h"
#include "fail.h"

#include <giac/global.h>
#include <giac/gausspol.h>

namespace GiNaC {

static unsigned int the_dimension = 7;

static giac::context * context_ptr=nullptr;
static giac::gen giac_zero;
static giac::gen giac_one;


inline giac::polynome gen2pol(const giac::gen& g) {
        return giac::polynome(giac::monomial<giac::gen>(g, the_dimension));
}

inline giac::gen num2gen(const numeric& n) {
        auto gp = n.to_giacgen(context_ptr);
        if (gp != nullptr)
                return std::move(*gp);
        else {
                std::stringstream ss;
                ss << n;
                return giac::gen(std::string(ss.str()), context_ptr);
        }
}

static giac::polynome replace_with_symbol(const ex& e, ex_int_map& map, exvector& revmap)
{
        // Expression already replaced? Then return the assigned symbol
        auto it = map.find(e);
        if (it != map.end()) {
                const int index = it->second;
                giac::monomial<giac::gen> mon(giac_one, index, the_dimension);
                return giac::polynome(mon);
        }

        // Otherwise create new symbol and add to dict
        const int index = revmap.size() + 1;
        map.insert(std::make_pair(e, index));
        revmap.push_back(e);
        giac::monomial<giac::gen> mon(giac_one, index, the_dimension);
        return giac::polynome(mon);
}

const giac::polynome basic::to_polynome(ex_int_map& map, exvector& revmap)
{
        throw std::runtime_error("basic::to_polynome: can't happen");
}

// Convert to giac polynomial over QQ, filling replacement dicts
// TODO: special case numeric mpz_t, int instead of string interface
const giac::polynome ex::to_polynome(ex_int_map& map, exvector& revmap) const
{
        if (is_exactly_a<add>(*this))
        {
                const add& a = ex_to<add>(*this);
                giac::polynome p = gen2pol(giac_zero);
                for (const auto& termex : a.seq) {
                        p = p + a.recombine_pair_to_ex(termex).to_polynome(map, revmap);
                }
                p = p + a.overall_coeff.to_polynome(map, revmap);
                return std::move(p);
        }
        else if (is_exactly_a<numeric>(*this))
        {
                const numeric& num = ex_to<numeric>(*this);
                if (num.is_real()) {
                        if (num.is_integer() or num.is_rational())
                                return gen2pol(num2gen(num));
                        else
                                return replace_with_symbol(num, map, revmap);
                } else { // complex
                        numeric re = num.real();
                        numeric im = num.imag();
                        giac::polynome re_p, im_p;
                        if (re.is_integer() or re.is_rational())
                                re_p = gen2pol(num2gen(re));
                        else
                                re_p = replace_with_symbol(re, map, revmap);
                        if (im.is_integer() or im.is_rational())
                                im_p = gen2pol(num2gen(im));
                        else
                                im_p = replace_with_symbol(im, map, revmap);
                        giac::polynome r = re_p + im_p * replace_with_symbol(I, map, revmap);
                        return std::move(r);
                }
        }
        else if (is_exactly_a<mul>(*this))
        {
                const mul& m = ex_to<mul>(*this);
                giac::polynome p = gen2pol(giac_one);
                for (const auto& termex : m.seq) {
                        p *= m.recombine_pair_to_ex(termex).to_polynome(map, revmap);
                }
                p *= m.overall_coeff.to_polynome(map, revmap);
                return std::move(p);
        }
        else if (is_exactly_a<power>(*this))
        {
                const power& pow = ex_to<power>(*this);
                if (is_exactly_a<numeric>(pow.exponent)) {
                        numeric expo = ex_to<numeric>(pow.exponent);
                        if (pow.exponent.info(info_flags::posint))
                                return std::move(giac::pow(pow.basis.to_polynome(map, revmap), expo.to_int()));
                        else if (pow.exponent.info(info_flags::negint))
                                return std::move(giac::pow(power(pow.basis, _ex_1).to_polynome(map, revmap), -expo.to_int()));
                }
                return replace_with_symbol(*this, map, revmap);
        }

        return replace_with_symbol(*this, map, revmap);
}

static ex gen2ex(const giac::gen& gen)
{
        // we need to handle giac types _INT_, _ZINT, and _CPLX
        switch (gen.type) {
                case giac::_INT_:
                        return numeric(gen.val);
                case giac::_ZINT:
                        mpz_t bigint;
                        mpz_init_set(bigint, *(gen.ref_ZINTptr()));
                        return numeric(bigint);
                case giac::_CPLX:
                        return gen2ex(gen.ref_CPLXptr()[0]) + I*gen2ex(gen.ref_CPLXptr()[1]);
                case giac::_FRAC:
                        return gen2ex(gen.ref_FRACptr()->num) /
                                                gen2ex(gen.ref_FRACptr()->den);
                default:
                        throw std::runtime_error("gen2ex: can't happen");
        }
}

static ex polynome_to_ex(const giac::polynome& p, const exvector& revmap)
{
        ex e = _ex0;
        for (const auto& mon : p.coord) {
                ex prod = _ex1;
                for (unsigned int varno=0; varno<mon.index.size(); ++varno) {
                        if (mon.index[varno] != 0)
                                prod *= power(revmap[varno], mon.index[varno]);
                }
                e += gen2ex(mon.value) * prod;
        }
        return e;
}

// GCD of two exes which are in polynomial form
ex gcdpoly(const ex &a, const ex &b, ex *ca=nullptr, ex *cb=nullptr, bool check_args=true)
{
//        std::cerr << "gcd(" << a << "," << b << ") = ";
        if (context_ptr == nullptr) {
                context_ptr=new giac::context();
                giac_zero = giac::gen(std::string("0"), context_ptr);
                giac_one = giac::gen(std::string("1"), context_ptr);
        }

        if (a.is_zero())
                return b;
        if (b.is_zero())
                return a;
        symbolset s1 = a.symbols();
        const symbolset& s2 = b.symbols();
        s1.insert(s2.begin(), s2.end());
        the_dimension = s1.size();

        ex_int_map map;
        exvector revmap;
        giac::polynome p = a.to_polynome(map, revmap);
        giac::polynome q = b.to_polynome(map, revmap);
        giac::polynome d(the_dimension);
        giac::gcd(p, q, d);
//        std::cerr << polynome_to_ex(d, revmap) << '\n';
        if (ca != nullptr) {
                giac::polynome quo;
                if (giac::exactquotient(p, d, quo))
                        *ca = polynome_to_ex(quo, revmap);
                else
                        throw(std::runtime_error("can't happen in gcdpoly"));
        }
        if (cb != nullptr) {
                giac::polynome quo;
                if (giac::exactquotient(q, d, quo))
                        *cb = polynome_to_ex(quo, revmap);
                else
                        throw(std::runtime_error("can't happen in gcdpoly"));
        }
        return polynome_to_ex(d, revmap);
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
		if (c.get_first_symbol(y))
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
	ex lcoef = r.coeff(x, deg);
	if (lcoef.info(info_flags::integer))
		return c;

	// GCD of all coefficients
	int ldeg = r.ldegree(x);
	if (deg == ldeg)
		return lcoef * c / lcoef.unit(x);
	ex cont = _ex0;
	for (int i=ldeg; i<=deg; i++)
		cont = gcdpoly(r.coeff(x, i), cont, nullptr, nullptr, false);
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
		if (is_exactly_a<symbol>(e.op(0)))
			return l;
		else {
			ex t = pow(lcmcoeff(e.op(0), l), ex_to<numeric>(e.op(1)));
                        if (is_exactly_a<numeric>(t))
                                return ex_to<numeric>(t);
                        else
                                return l;
                }
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
numeric lcm_of_coefficients_denominators(const ex &e)
{
	return lcmcoeff(e, *_num1_p);
}

/** Bring polynomial from Q[X] to Z[X] by multiplying in the previously
 *  determined LCM of the coefficient's denominators.
 *
 *  @param e  multivariate polynomial (need not be expanded)
 *  @param lcm  LCM to multiply in */
ex multiply_lcm(const ex &e, const numeric &lcm)
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
		if (is_exactly_a<symbol>(e.op(0)))
			return e * lcm;
		else {
			ex root_of_lcm = lcm.power(ex_to<numeric>(e.op(1)).inverse());
			if (is_exactly_a<numeric>(root_of_lcm)
                                        and ex_to<numeric>(root_of_lcm).is_rational())
				return pow(multiply_lcm(e.op(0), ex_to<numeric>(root_of_lcm)), e.op(1));
			else
				return e * lcm;
		}
	} else
		return e * lcm;
}


} // namespace GiNaC

#endif // HAVE_LIBGIAC
