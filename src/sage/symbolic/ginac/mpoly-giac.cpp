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

} // namespace GiNaC

#endif // HAVE_LIBGIAC
