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
        numeric numn = ex_to<numeric>(n);
        numeric numx = ex_to<numeric>(x);
        std::vector<numeric> numveca, numvecb;
        numveca.push_back(numn / *_num_2_p);
        numveca.push_back(*_num1_2_p + (numn / *_num_2_p));
        return hypergeometric_2F1(numveca, numvecb, -pow(numx, *_num2_p).inverse(), parent);
}

static ex hermite_eval(const ex& n, const ex& x)
{
	if (is_exactly_a<numeric>(x)) {
                numeric numx = ex_to<numeric>(x);
                if (numx.is_zero())
                        if (n.info(info_flags::integer)
                                and n.info(info_flags::odd))
                                return _ex0;
                if (is_exactly_a<numeric>(n) and not numx.info(info_flags::crational))
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
        numeric coeff = _num2_p->power(numn);
        ex sum = _ex0;
        int fac = 1;
        while (numn >= 0) {
                sum = sum + power(x, numn) * coeff;
                coeff /= *_num_4_p;
                coeff *= numn;
                --numn;
                coeff *= numn;
                --numn;
                coeff /= fac++;
                }
        return sum;
}

static ex hermite_deriv(const ex& n, const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	return _ex2 * n * hermite(n-1, x).hold();
}

REGISTER_FUNCTION(hermite, eval_func(hermite_eval).
                        evalf_func(hermite_evalf).
                        derivative_func(hermite_deriv).
		        latex_name("H"));

} // namespace GiNaC
