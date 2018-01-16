/** @file inifcns_orthopoly.cpp
 *
 *  Implementation of orthogonal polynomials as functions.
 *
 *  (C) 2016 Ralf Stephan <ralf@ark.in-berlin.de>
 */

#include "inifcns.h"
#include "ex.h"
#include "constant.h"
#include "infinity.h"
#include "numeric.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "symbol.h"
#include "pseries.h"
#include "utils.h"

namespace GiNaC {

//////////
// Hermite polynomials H_n(x)
//////////

static ex hermite_evalf(const ex& n, const ex& x, PyObject* parent)
{
	if (not is_exactly_a<numeric>(x)
             or not is_exactly_a<numeric>(n))
                return hermite(n,x).hold();

        // see http://dlmf.nist.gov/18.5.E13
        numeric numn = ex_to<numeric>(n);
        numeric numx = ex_to<numeric>(x);
        std::vector<numeric> numveca, numvecb;
        numveca.push_back(numn / *_num_2_p);
        numveca.push_back(*_num1_2_p + (numn / *_num_2_p));
        return pow(numx * (*_num2_p), numn) * hypergeometric_pFq(numveca, numvecb, numx.power(-2).negative(), parent);
}

static ex hermite_eval(const ex& n, const ex& x)
{
	if (is_exactly_a<numeric>(x)) {
                numeric numx = ex_to<numeric>(x);
                if (numx.is_zero())
                        if (n.info(info_flags::integer)
                                and n.info(info_flags::odd))
                                return _ex0;
                if (is_exactly_a<numeric>(n) and numx.info(info_flags::inexact))
                        return hermite_evalf(n, x, nullptr);
        }

        if (not is_exactly_a<numeric>(n))
                return hermite(n, x).hold();
        numeric numn = ex_to<numeric>(n);
        if (not numn.info(info_flags::integer) or numn < 0)
                throw std::runtime_error("hermite_eval: The index n must be a nonnegative integer");

        if (numn.is_zero())
                return _ex1;

        // apply the formula from http://oeis.org/A060821
        // T(n, k) = ((-1)^((n-k)/2))*(2^k)*n!/(k!*((n-k)/2)!) if n-k is even and >=0, else 0.
        // sequentially for all viable k. Effectively there is the recurrence
        // T(n, k) = -(k+2)*(k+1)/(2*(n-k)) * T(n, k+2), with T(n, n) = 2^n
        numeric coeff = _num2_p->pow_intexp(numn);
        ex sum = _ex0;
        int fac = 1;
        epvector vec;
        while (numn >= 0) {
                vec.push_back(expair(power(x, numn), coeff));
                coeff /= *_num_4_p;
                coeff *= numn;
                --numn;
                coeff *= numn;
                --numn;
                coeff /= fac++;
                }
        return add(vec);
}

static ex hermite_deriv(const ex& n, const ex & x, unsigned deriv_param)
{
	    if (deriv_param == 0)
                    throw std::runtime_error("derivative w.r.t. to the index is not supported yet");
	    return _ex2 * n * hermite(n-1, x).hold();
}

REGISTER_FUNCTION(hermite, eval_func(hermite_eval).
                        evalf_func(hermite_evalf).
                        derivative_func(hermite_deriv).
		        latex_name("H"));

//////////
// Gegenbauer (ultraspherical) polynomials C^a_n(x)
//////////

static ex gegenb_evalf(const ex& n, const ex &a, const ex& x, PyObject* parent)
{
	if (not is_exactly_a<numeric>(x)
             or not is_exactly_a<numeric>(a)
             or not is_exactly_a<numeric>(n))
                return gegenbauer(n,a,x).hold();

        // see http://dlmf.nist.gov/18.5.E9
        numeric numn = ex_to<numeric>(n);
        numeric numx = ex_to<numeric>(x);
        numeric numa = ex_to<numeric>(a);
        numeric num2a = numa * (*_num2_p);
        std::vector<numeric> numveca, numvecb;
        numveca.push_back(-numn);
        numveca.push_back(numn + num2a);
        numvecb.push_back(numa + (*_num1_2_p));
        numeric factor = numn + num2a - (*_num1_p);
        factor = factor.binomial(numn);
        return factor * hypergeometric_pFq(numveca, numvecb, ((*_num1_p)-numx)/(*_num2_p), parent);
}

static ex gegenb_eval(const ex& n, const ex &a, const ex& x)
{
	numeric numa = 1;

	if (is_exactly_a<numeric>(x)) {
                numeric numx = ex_to<numeric>(x);
                if (is_exactly_a<numeric>(n)
			and is_exactly_a<numeric>(a)
			and numx.info(info_flags::inexact))
                        return gegenb_evalf(n, a, x, nullptr);
        }

        if (not is_exactly_a<numeric>(n))
                return gegenbauer(n, a, x).hold();

        numeric numn = ex_to<numeric>(n);
        if (not numn.info(info_flags::integer) or numn < 0)
                throw std::runtime_error("gegenb_eval: The index n must be a nonnegative integer");

        if (numn.is_zero())
                return _ex1;
	if (numn.is_equal(*_num1_p))
		return _ex2 * a * x;

	if (is_exactly_a<numeric>(a)) {
		numa = ex_to<numeric>(a);
		if (numa.is_zero())
			return _ex0; // C_n^0(x) = 0
	}

	if (not is_exactly_a<numeric>(a) or numa < 0) {
		ex p = _ex1;
		ex sum = _ex0;
		ex sign = _ex1;
		ex aa = _ex1;
		long nn = numn.to_long();

		sign = - 1;
		for (long k=0; k<=nn/2; k++) {
			sign *= - 1;
			aa = a;
			p = 1;

			// rising factorial (a)_{n-k}
			for (long i=0; i<nn-k; i++) {
				p *= aa;
				aa = aa + 1;
			}

			p /= numeric(k).factorial();
			p /= numeric(nn-2*k).factorial();
			sum += pow(2*x, nn-2*k) * p * sign;
		}
		return sum;
	}

	// from here on see flint2/fmpq_poly/gegenbauer_c.c
	numeric numer = numa.numer();
	numeric denom = numa.denom();
	numeric t = numn.factorial();
	numeric overall_denom = denom.pow_intexp(numn) * t;

	long nn = numn.to_long();
	numeric p = t / (numeric(nn/2).factorial());
	if ((nn%2) != 0u)
		p *= *_num2_p;
	if ((nn&2) != 0u)
		p = -p;

	for (long k=0; k < nn-nn/2; k++) {
		p *= numer;
		numer += denom;
	}

	p *= denom.power(nn/2);
        ex sum = _ex0;
	if ((nn%2) != 0u)
		sum += x*p;
	else
		sum += p;

	for (long k=nn/2-1; k>=0; --k) {
		p *= numer;
		p *= 4*(k+1);
		p *= *_num_1_p;
		p /= denom;
		p /= nn-2*k-1;
		p /= nn-2*k;
		sum += pow(x, nn-2*k) * p;
		numer += denom;
	}

        return sum / overall_denom;
}

static ex gegenb_deriv(const ex& n, const ex & a, const ex & x, unsigned deriv_param)
{
	    if (deriv_param == 0)
                    throw std::runtime_error("derivative w.r.t. to the index is not supported yet");
	    if (deriv_param == 1)
                    throw std::runtime_error("derivative w.r.t. to the second index is not supported yet");

	    return _ex2 * a * gegenbauer(n-1, a+1, x).hold();
}

REGISTER_FUNCTION(gegenbauer, eval_func(gegenb_eval).
                        evalf_func(gegenb_evalf).
                        derivative_func(gegenb_deriv).
		        latex_name("C"));

} // namespace GiNaC
