/** @file mpoly-singular.cpp
 *
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
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

#include <factory/factory.h>

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
#include "relational.h"
#include "symbol.h"
#include "function.h"
#include "utils.h"
#include "fail.h"

namespace GiNaC {


static CanonicalForm num2canonical(const numeric& n)
{
        if (n.is_mpz() or n.is_mpq())
                return n.to_canonical();

        std::cerr << "\nn = " << n << "\n";
        std::cerr << "is_pyobject = " << n.is_pyobject() << "\n";
        throw std::runtime_error("num2canonical: can't happen");
}

static CanonicalForm replace_with_symbol(const ex& e, ex_int_map& map, exvector& revmap)
{
        // Expression already replaced? Then return the assigned symbol
        auto it = map.find(e);
        if (it != map.end()) {
//std::cerr<< "found Var("<<it->second<<")\n";
                return Variable(it->second);
        }

        // Otherwise create new symbol and add to dict
        const int index = revmap.size() + 1;
        map.insert(std::make_pair(e, index));
//std::cerr<<"\n"<<e<<" not found; added as Var("<<index<<")\n";
        revmap.push_back(e);
        return Variable(index);
}

const CanonicalForm basic::to_canonical(ex_int_map& map, exvector& revmap)
{
        throw std::runtime_error("basic::to_polynome: can't happen");
}

// Convert to Singulat polynomial over ZZ, filling replacement dicts
// TODO: special case numeric mpz_t, int instead of string interface
const CanonicalForm ex::to_canonical(ex_int_map& map, exvector& revmap) const
{
        if (is_exactly_a<add>(*this))
        {
                const add& a = ex_to<add>(*this);
                CanonicalForm p(0);
                for (const auto& termex : a.seq) {
                        p = p + a.recombine_pair_to_ex(termex).to_canonical(map, revmap);
                }
                p = p + a.overall_coeff.to_canonical(map, revmap);
                return p;
        }
        else if (is_exactly_a<numeric>(*this))
        {
                const numeric& num = ex_to<numeric>(*this);
                if (num.is_real()) {
                        if (num.is_integer() or num.is_rational())
                                return num2canonical(num);
                        else
                                return replace_with_symbol(num, map, revmap);
                } else { // complex
                        numeric re = num.real();
                        numeric im = num.imag();
                        CanonicalForm re_p, im_p;
                        if (re.is_integer() or re.is_rational())
                                re_p = num2canonical(re);
                        else
                                re_p = replace_with_symbol(re, map, revmap);
                        if (im.is_integer() or im.is_rational())
                                im_p = num2canonical(im);
                        else
                                im_p = replace_with_symbol(im, map, revmap);
                        return re_p + im_p * replace_with_symbol(I, map, revmap);
                }
        }
        else if (is_exactly_a<mul>(*this))
        {
                const mul& m = ex_to<mul>(*this);
                CanonicalForm p(1);
                for (const auto& termex : m.seq) {
                        p = p * m.recombine_pair_to_ex(termex).to_canonical(map, revmap);
                }
                p *= m.overall_coeff.to_canonical(map, revmap);
                return p;
        }
        else if (is_exactly_a<power>(*this))
        {
                const power& pow = ex_to<power>(*this);
                if (is_exactly_a<numeric>(pow.exponent)) {
                        numeric expo = ex_to<numeric>(pow.exponent);
                        if (pow.exponent.info(info_flags::posint))
                                return ::power(pow.basis.to_canonical(map, revmap), expo.to_int());
                        else if (pow.exponent.info(info_flags::negint))
                                return ::power(GiNaC::power(pow.basis, _ex_1).to_canonical(map, revmap), -expo.to_int());
                }
                return replace_with_symbol(*this, map, revmap);
        }

        return replace_with_symbol(*this, map, revmap);
}

static numeric can2num(const CanonicalForm& f)
{
        if (f.isOne())
                return *_num1_p;
        if (f.isImm())
                return numeric(f.intval());
        mpz_t bignum;
        f.mpzval(bignum);
        return numeric(bignum);
}

static ex coeff_to_ex(const CanonicalForm& f, const exvector& revmap)
{
        if (f.isOne())
                return _ex1;
        if (f.isImm())
                return numeric(f.intval());
        if (f.inZ())
                throw std::runtime_error("can't happen #1");
        if (f.inQ()) {
                CanonicalForm num = f.num();
                CanonicalForm den = f.den();
                return can2num(num)/can2num(den);
        }
        ex res = _ex0;
        for ( CFIterator I = f; I.hasTerms(); I++ ) {
//std::cerr<<"Lev: "<<f.level()<<", exp: "<<i.exp()<<"\n";
                res += mul(GiNaC::power(revmap[f.level()-1], I.exp()),
                                 coeff_to_ex(I.coeff(), revmap));
        }
        return res;
}

static ex canonical_to_ex(const CanonicalForm& f, const exvector& revmap)
{
        if (f.isOne())
                return _ex1;
        if (f.inCoeffDomain()) {
                if (f.isImm())
                        return numeric(f.intval());
                if (f.inZ())
                        throw std::runtime_error("can't happen #2");
                        mpz_t bigint;
                        f.mpzval(bigint);
                        numeric n(bigint);
                        mpz_clear(bigint);
                        return n;
                if (f.inQ()) {
                        CanonicalForm num = f.num();
                        CanonicalForm den = f.den();
                        mpz_t bigintnum, bigintden;
                        num.mpzval(bigintnum);
                        den.mpzval(bigintden);
                        numeric n(bigintnum);
                        numeric d(bigintden);
                        mpz_clear(bigintnum);
                        mpz_clear(bigintden);
                        return n/d;
                }
                else {
                        throw std::runtime_error("can't happen #3");
                }
        }

        ex e = _ex0;
        for ( CFIterator i = f; i.hasTerms(); i++ ) {
//std::cerr<<"Lev: "<<f.level()<<", exp: "<<i.exp()<<"\n";
                e += mul(coeff_to_ex(i.coeff(), revmap),
                                   power(revmap[f.level()-1], i.exp()));
        }
        return e;
}

// GCD of two exes which are in polynomial form
// If giac is requested we stand back
#ifndef PYNAC_HAVE_LIBGIAC
ex gcdpoly(const ex &a, const ex &b, ex *ca=nullptr, ex *cb=nullptr, bool check_args=true)
{
//std::cerr << "gcd(" << a << "," << b << ") = ";

        if (a.is_zero())
                return b;
        if (b.is_zero())
                return a;

        // GCD of numerics
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b)) {
		numeric g = gcd(ex_to<numeric>(a), ex_to<numeric>(b));
		if ((ca != nullptr) || (cb != nullptr)) {
			if (g.is_zero()) {
				if (ca != nullptr)
					*ca = _ex0;
				if (cb != nullptr)
					*cb = _ex0;
			} else {
				if (ca != nullptr)
					*ca = ex_to<numeric>(a) / g;
				if (cb != nullptr)
					*cb = ex_to<numeric>(b) / g;
			}
		}
		return g;
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
			g.push_back(gcdpoly(a.op(i), part_b, &part_ca, &part_cb, check_args));
			acc_ca.push_back(part_ca);
			part_b = part_cb;
		}
		if (ca != nullptr)
			*ca = (new mul(acc_ca))->setflag(status_flags::dynallocated);
		if (cb != nullptr)
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
			g.push_back(gcdpoly(part_a, b.op(i), &part_ca, &part_cb, check_args));
			acc_cb.push_back(part_cb);
			part_a = part_ca;
		}
		if (ca != nullptr)
			*ca = part_a;
		if (cb != nullptr)
			*cb = (new mul(acc_cb))->setflag(status_flags::dynallocated);
		return (new mul(g))->setflag(status_flags::dynallocated);
	}

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
					if (ca != nullptr)
						*ca = _ex1;
					if (cb != nullptr)
						*cb = power(p, exp_b - exp_a);
					return power(p, exp_a);
				} else {
					if (ca != nullptr)
						*ca = power(p, exp_a - exp_b);
					if (cb != nullptr)
						*cb = _ex1;
					return power(p, exp_b);
				}
			} else {
				ex p_co, pb_co;
				ex p_gcd = gcdpoly(p, pb, &p_co, &pb_co, check_args);
				if (p_gcd.is_equal(_ex1)) {
					// a(x) = p(x)^n, b(x) = p_b(x)^m, gcd (p, p_b) = 1 ==>
					// gcd(a,b) = 1
					if (ca != nullptr)
						*ca = a;
					if (cb != nullptr)
						*cb = b;
					return _ex1;
					// XXX: do I need to check for p_gcd = -1?
				} else {
					// there are common factors:
					// a(x) = g(x)^n A(x)^n, b(x) = g(x)^m B(x)^m ==>
					// gcd(a, b) = g(x)^n gcd(A(x)^n, g(x)^(n-m) B(x)^m
					if (exp_a < exp_b) {
						return power(p_gcd, exp_a)*
							gcdpoly(power(p_co, exp_a), power(p_gcd, exp_b-exp_a)*power(pb_co, exp_b), ca, cb, false);
					} else {
						return power(p_gcd, exp_b)*
							gcdpoly(power(p_gcd, exp_a - exp_b)*power(p_co, exp_a), power(pb_co, exp_b), ca, cb, false);
					}
				} // p_gcd.is_equal(_ex1)
			} // p.is_equal(pb)

		} else {
			if (p.is_equal(b)) {
				// a = p^n, b = p, gcd = p
				if (ca != nullptr)
					*ca = power(p, a.op(1) - 1);
				if (cb != nullptr)
					*cb = _ex1;
				return p;
			} 

			ex p_co, bpart_co;
			ex p_gcd = gcdpoly(p, b, &p_co, &bpart_co, false);

			if (p_gcd.is_equal(_ex1)) {
				// a(x) = p(x)^n, gcd(p, b) = 1 ==> gcd(a, b) = 1
				if (ca != nullptr)
					*ca = a;
				if (cb != nullptr)
					*cb = b;
				return _ex1;
			} else {
				// a(x) = g(x)^n A(x)^n, b(x) = g(x) B(x) ==> gcd(a, b) = g(x) gcd(g(x)^(n-1) A(x)^n, B(x))
				return p_gcd*gcdpoly(power(p_gcd, exp_a-1)*power(p_co, exp_a), bpart_co, ca, cb, false);
			}
		} // is_exactly_a<power>(b)

	} else if (is_exactly_a<power>(b)) {
		ex p = b.op(0);
		if (p.is_equal(a)) {
			// a = p, b = p^n, gcd = p
			if (ca != nullptr)
				*ca = _ex1;
			if (cb != nullptr)
				*cb = power(p, b.op(1) - 1);
			return p;
		}

		ex p_co, apart_co;
		const ex& exp_b(b.op(1));
		ex p_gcd = gcdpoly(a, p, &apart_co, &p_co, false);
		if (p_gcd.is_equal(_ex1)) {
			// b=p(x)^n, gcd(a, p) = 1 ==> gcd(a, b) == 1
			if (ca != nullptr)
				*ca = a;
			if (cb != nullptr)
				*cb = b;
			return _ex1;
		} else {
			// there are common factors:
			// a(x) = g(x) A(x), b(x) = g(x)^n B(x)^n ==> gcd = g(x) gcd(g(x)^(n-1) A(x)^n, B(x))

			return p_gcd*gcdpoly(apart_co, power(p_gcd, exp_b-1)*power(p_co, exp_b), ca, cb, false);
		} // p_gcd.is_equal(_ex1)
	}


        ex_int_map map;
        exvector revmap;
        On(SW_RATIONAL);
        setCharacteristic(0);
        CanonicalForm p = a.to_canonical(map, revmap);
        CanonicalForm q = b.to_canonical(map, revmap);
        CanonicalForm d = gcd(p, q);

//std::cerr << canonical_to_ex(d, revmap) << '\n';
        if (ca != nullptr) {
                CanonicalForm quo;
                if (fdivides(d, p, quo))
                        *ca = canonical_to_ex(quo, revmap);
                else
                        throw(std::runtime_error("can't happen in gcdpoly"));
        }
        if (cb != nullptr) {
                CanonicalForm quo;
                if (fdivides(d, q, quo))
                        *cb = canonical_to_ex(quo, revmap);
                else
                        throw(std::runtime_error("can't happen in gcdpoly"));
        }
        ex res = canonical_to_ex(d, revmap);
        return res;
}
#endif //PYNAC_HAVE_LIBGIAC

} // namespace GiNaC

