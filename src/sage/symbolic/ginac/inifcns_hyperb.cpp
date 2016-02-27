/** @file inifcns_hyperb.cpp
 *
 *  Implementation of hyperbolic functions.
 *
 *  (C) 2016 Ralf Stephan <ralf@ark.in-berlin.de>
 */

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

#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <memory>

namespace GiNaC {


/* In Sage all the area hyperbolic functions are printed with "arc" instead
   of "a" at the beginning.   This is for consistency with other
   computer algebra systems.   These print methods are registered
   below with each of the corresponding inverse trig function. */


//////////
// hyperbolic sine (trigonometric function)
//////////

static ex sinh_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return sinh(ex_to<numeric>(x));
	
	return sinh(x).hold();
}

static ex sinh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// sinh(0) -> 0
		if (x.is_zero())
			return _ex0;        

		// sinh(float) -> float
		if (!x.info(info_flags::crational))
			return sinh(ex_to<numeric>(x));

		// sinh() is odd
		if (x.info(info_flags::negative))
			return -sinh(-x);
	}

	// sinh(oo) -> oo
	// sinh(-oo) -> -oo
	// sinh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(UnsignedInfinity))
			throw (std::runtime_error("sinh_eval(): sinh(unsigned_infinity) encountered"));
		return x;
	}
	
        ex xoverpi = x/Pi;
	if (is_exactly_a<numeric>(xoverpi) &&
		ex_to<numeric>(xoverpi).real().is_zero())  // sinh(I*x) -> I*sin(x)
		return I*sin(x/I);
	
	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// sinh(asinh(x)) -> x
		if (is_ex_the_function(x, asinh))
			return t;

		// sinh(acosh(x)) -> sqrt(x-1) * sqrt(x+1)
		if (is_ex_the_function(x, acosh))
			return sqrt(t-_ex1)*sqrt(t+_ex1);

		// sinh(atanh(x)) -> x/sqrt(1-x^2)
		if (is_ex_the_function(x, atanh))
			return t*power(_ex1-power(t,_ex2),_ex_1_2);
	}
	
	return sinh(x).hold();
}

static ex sinh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx sinh(x) -> cosh(x)
	return cosh(x);
}

static ex sinh_real_part(const ex & x)
{
	return sinh(GiNaC::real_part(x))*cos(GiNaC::imag_part(x));
}

static ex sinh_imag_part(const ex & x)
{
	return cosh(GiNaC::real_part(x))*sin(GiNaC::imag_part(x));
}

static ex sinh_conjugate(const ex & x)
{
	// conjugate(sinh(x))==sinh(conjugate(x))
	return sinh(x.conjugate());
}

REGISTER_FUNCTION(sinh, eval_func(sinh_eval).
                        evalf_func(sinh_evalf).
                        derivative_func(sinh_deriv).
                        real_part_func(sinh_real_part).
                        imag_part_func(sinh_imag_part).
			conjugate_func(sinh_conjugate).
		        latex_name("\\sinh"));

//////////
// hyperbolic cosine (trigonometric function)
//////////

static ex cosh_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return cosh(ex_to<numeric>(x));
	
	return cosh(x).hold();
}

static ex cosh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// cosh(0) -> 1
		if (x.is_zero())
			return _ex1;

		// cosh(float) -> float
		if (!x.info(info_flags::crational))
			return cosh(ex_to<numeric>(x));

		// cosh() is even
		if (x.info(info_flags::negative))
			return cosh(-x);
	}
	
	// cosh(oo) -> oo
	// cosh(-oo) -> oo
	// cosh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(UnsignedInfinity))
			throw (std::runtime_error("cosh_eval(): cosh(unsigned_infinity) encountered"));
		return Infinity;
	}

        ex xoverpi = x/Pi;
	if (is_exactly_a<numeric>(xoverpi) &&
		ex_to<numeric>(xoverpi).real().is_zero())  // cosh(I*x) -> cos(x)
		return cos(x/I);
	
	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// cosh(acosh(x)) -> x
		if (is_ex_the_function(x, acosh))
			return t;

		// cosh(asinh(x)) -> sqrt(1+x^2)
		if (is_ex_the_function(x, asinh))
			return sqrt(_ex1+power(t,_ex2));

		// cosh(atanh(x)) -> 1/sqrt(1-x^2)
		if (is_ex_the_function(x, atanh))
			return power(_ex1-power(t,_ex2),_ex_1_2);
	}
	
	return cosh(x).hold();
}

static ex cosh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx cosh(x) -> sinh(x)
	return sinh(x);
}

static ex cosh_real_part(const ex & x)
{
	return cosh(GiNaC::real_part(x))*cos(GiNaC::imag_part(x));
}

static ex cosh_imag_part(const ex & x)
{
	return sinh(GiNaC::real_part(x))*sin(GiNaC::imag_part(x));
}

static ex cosh_conjugate(const ex & x)
{
	// conjugate(cosh(x))==cosh(conjugate(x))
	return cosh(x.conjugate());
}

REGISTER_FUNCTION(cosh, eval_func(cosh_eval).
                        evalf_func(cosh_evalf).
                        derivative_func(cosh_deriv).
                        real_part_func(cosh_real_part).
                        imag_part_func(cosh_imag_part).
                        conjugate_func(cosh_conjugate).
                        latex_name("\\cosh"));

//////////
// hyperbolic tangent (trigonometric function)
//////////

static ex tanh_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return tanh(ex_to<numeric>(x));
	
	return tanh(x).hold();
}

static ex tanh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// tanh(0) -> 0
		if (x.is_zero())
			return _ex0;

		// tanh(float) -> float
		if (!x.info(info_flags::crational))
			return tanh(ex_to<numeric>(x));

		// tanh() is odd
		if (x.info(info_flags::negative))
			return -tanh(-x);
	}
	
	// tanh(oo) -> 1
	// tanh(-oo) -> -1
	// tanh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity))
			return _ex1;
		if (x.is_equal(NegInfinity))
			return _ex_1;
		// x is UnsignedInfinity
		throw (std::runtime_error("tanh_eval(): tanh(unsigned_infinity) encountered"));
	}
		
        ex xoverpi = x/Pi;
	if (is_exactly_a<numeric>(xoverpi) &&
		ex_to<numeric>(xoverpi).real().is_zero())  // tanh(I*x) -> I*tan(x);
		return I*tan(x/I);
	
	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// tanh(atanh(x)) -> x
		if (is_ex_the_function(x, atanh))
			return t;

		// tanh(asinh(x)) -> x/sqrt(1+x^2)
		if (is_ex_the_function(x, asinh))
			return t*power(_ex1+power(t,_ex2),_ex_1_2);

		// tanh(acosh(x)) -> sqrt(x-1)*sqrt(x+1)/x
		if (is_ex_the_function(x, acosh))
			return sqrt(t-_ex1)*sqrt(t+_ex1)*power(t,_ex_1);
	}
	
	return tanh(x).hold();
}

static ex tanh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx tanh(x) -> 1-tanh(x)^2
	return _ex1-power(tanh(x),_ex2);
}

static ex tanh_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole falls back to tanh_deriv.
	// On a pole simply expand sinh(x)/cosh(x).
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (!(2*I*x_pt/Pi).info(info_flags::odd))
		throw do_taylor();  // caught by function::series()
	// if we got here we have to care for a simple pole
	return (sinh(x)/cosh(x)).series(rel, order, options);
}

// See http://dlmf.nist.gov/4.35.E36
static ex tanh_real_part(const ex & x)
{
	ex a = mul(real_part(x), _ex2);
	ex b = mul(imag_part(x), _ex2);
	return sinh(a) / (cosh(a) + cos(b));
}

static ex tanh_imag_part(const ex & x)
{
	ex a = mul(real_part(x), _ex2);
	ex b = mul(imag_part(x), _ex2);
	return sin(b) / (cosh(a) + cos(b));
}

static ex tanh_conjugate(const ex & x)
{
	// conjugate(tanh(x))==tanh(conjugate(x))
	return tanh(x.conjugate());
}

REGISTER_FUNCTION(tanh, eval_func(tanh_eval).
                        evalf_func(tanh_evalf).
                        derivative_func(tanh_deriv).
                        series_func(tanh_series).
                        real_part_func(tanh_real_part).
                        imag_part_func(tanh_imag_part).
                        conjugate_func(tanh_conjugate).
                        latex_name("\\tanh"));

//////////
// inverse hyperbolic sine (trigonometric function)
//////////

static ex asinh_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return asinh(ex_to<numeric>(x));
	
	return asinh(x).hold();
}

static ex asinh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// asinh(0) -> 0
		if (x.is_zero())
			return _ex0;

		// asinh(float) -> float
		if (!x.info(info_flags::crational))
			return asinh(ex_to<numeric>(x));

		// asinh() is odd
		if (x.info(info_flags::negative))
			return -asinh(-x);
	}
	
	// asinh(oo) -> oo
	// asinh(-oo) -> -oo
	// asinh(UnsignedInfinity) -> UnsignedInfinity
	if (x.info(info_flags::infinity)) {
		return x;
	}

	return asinh(x).hold();
}

static ex asinh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx asinh(x) -> 1/sqrt(1+x^2)
	return power(_ex1+power(x,_ex2),_ex_1_2);
}

static ex asinh_conjugate(const ex & x)
{
	// conjugate(asinh(x))==asinh(conjugate(x)) unless on the branch cuts which
	// run along the imaginary axis outside the interval [-I, +I].
	if (x.info(info_flags::real))
		return asinh(x);
	if (is_exactly_a<numeric>(x)) {
		const numeric x_re = ex_to<numeric>(x.real_part());
		const numeric x_im = ex_to<numeric>(x.imag_part());
		if (!x_re.is_zero() ||
		    (x_im > *_num_1_p && x_im < *_num1_p))
			return asinh(x.conjugate());
	}
	return conjugate_function(asinh(x)).hold();
}

REGISTER_FUNCTION(asinh, eval_func(asinh_eval).
                         evalf_func(asinh_evalf).
                         derivative_func(asinh_deriv).
                         conjugate_func(asinh_conjugate).
			 set_name("arcsinh"));

//////////
// inverse hyperbolic cosine (trigonometric function)
//////////

static ex acosh_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return acosh(ex_to<numeric>(x));
	
	return acosh(x).hold();
}

static ex acosh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// acosh(0) -> Pi*I/2
		if (x.is_zero())
			return Pi*I*numeric(1,2);

		// acosh(1) -> 0
		if (x.is_equal(_ex1))
			return _ex0;

		// acosh(-1) -> Pi*I
		if (x.is_equal(_ex_1))
			return Pi*I;

		// acosh(float) -> float
		if (!x.info(info_flags::crational))
			return acosh(ex_to<numeric>(x));

		// acosh(-x) -> Pi*I-acosh(x)
		if (x.info(info_flags::negative))
			return Pi*I-acosh(-x);
	}
	
	// acosh(oo) -> oo
	// acosh(-oo) -> oo
	// acosh(UnsignedInfinity) -> oo
	if (x.info(info_flags::infinity)) {
		return Infinity;
	}

	return acosh(x).hold();
}

static ex acosh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx acosh(x) -> 1/(sqrt(x-1)*sqrt(x+1))
	return power(x+_ex_1,_ex_1_2)*power(x+_ex1,_ex_1_2);
}

static ex acosh_conjugate(const ex & x)
{
	// conjugate(acosh(x))==acosh(conjugate(x)) unless on the branch cut
	// which runs along the real axis from +1 to -inf.
	if (is_exactly_a<numeric>(x) &&
	    (!x.imag_part().is_zero() || x > *_num1_p)) {
		return acosh(x.conjugate());
	}
	return conjugate_function(acosh(x)).hold();
}

REGISTER_FUNCTION(acosh, eval_func(acosh_eval).
                         evalf_func(acosh_evalf).
                         derivative_func(acosh_deriv).
                         conjugate_func(acosh_conjugate).
			 set_name("arccosh"));

//////////
// inverse hyperbolic tangent (trigonometric function)
//////////

static ex atanh_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return atanh(ex_to<numeric>(x));
	
	return atanh(x).hold();
}

static ex atanh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// atanh(0) -> 0
		if (x.is_zero())
			return _ex0;

		/*
		// atanh({+|-}1) -> throw
		if (x.is_equal(_ex1) || x.is_equal(_ex_1))
			throw (pole_error("atanh_eval(): logarithmic pole",0));
		*/

		// atanh(1) -> oo
		if (x.is_equal(_ex1))
			return Infinity;
		// atahn(-1) -> -oo
		if (x.is_equal(_ex_1))
			return NegInfinity;

		// atanh(float) -> float
		if (!x.info(info_flags::crational))
			return atanh(ex_to<numeric>(x));

		// atanh() is odd
		if (x.info(info_flags::negative))
			return -atanh(-x);
	}
	
	// atanh(oo) -> -i*pi/2
	// atanh(-oo) -> i*pi/2
	// atanh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity))
			return _ex_1_2*Pi*I;
		if (x.is_equal(NegInfinity))
			return _ex1_2*Pi*I;
		// x is UnsignedInfinity
		throw (std::runtime_error("arctanh_eval(): arctanh(unsigned_infinity) encountered"));
	}
		
	return atanh(x).hold();
}

static ex atanh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx atanh(x) -> 1/(1-x^2)
	return power(_ex1-power(x,_ex2),_ex_1);
}

static ex atanh_series(const ex &arg,
                       const relational &rel,
                       int order,
                       unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole or cut falls back to atanh_deriv.
	// There are two branch cuts, one runnig from 1 up the real axis and one
	// one running from -1 down the real axis.  The points 1 and -1 are poles
	// On the branch cuts and the poles series expand
	//     (log(1+x)-log(1-x))/2
	// instead.
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (!(arg_pt).info(info_flags::real))
		throw do_taylor();     // Im(x) != 0
	if ((arg_pt).info(info_flags::real) && abs(arg_pt)<_ex1)
		throw do_taylor();     // Im(x) == 0, but abs(x)<1
	// care for the poles, using the defining formula for atanh()...
	if (arg_pt.is_equal(_ex1) || arg_pt.is_equal(_ex_1))
		return ((log(_ex1+arg)-log(_ex1-arg))*_ex1_2).series(rel, order, options);
	// ...and the branch cuts (the discontinuity at the cut being just I*Pi)
	if (!(options & series_options::suppress_branchcut)) {
 		// method:
 		// This is the branch cut: assemble the primitive series manually and
 		// then add the corresponding complex step function.
 		const symbol &s = ex_to<symbol>(rel.lhs());
 		const ex &point = rel.rhs();
 		const symbol foo;
 		const ex replarg = series(atanh(arg), s==foo, order).subs(foo==point, subs_options::no_pattern);
		ex Order0correction = replarg.op(0)+csgn(I*arg)*Pi*I*_ex1_2;
		if (arg_pt<_ex0)
			Order0correction += log((arg_pt+_ex_1)/(arg_pt+_ex1))*_ex1_2;
		else
			Order0correction += log((arg_pt+_ex1)/(arg_pt+_ex_1))*_ex_1_2;
 		epvector seq;
		seq.push_back(expair(Order0correction, _ex0));
 		seq.push_back(expair(Order(_ex1), order));
 		return series(replarg - pseries(rel, seq), rel, order);
	}
	throw do_taylor();
}

static ex atanh_conjugate(const ex & x)
{
	// conjugate(atanh(x))==atanh(conjugate(x)) unless on the branch cuts which
	// run along the real axis outside the interval [-1, +1].
	if (is_exactly_a<numeric>(x) &&
	    (!x.imag_part().is_zero() || (x > *_num_1_p && x < *_num1_p))) {
		return atanh(x.conjugate());
	}
	return conjugate_function(atanh(x)).hold();
}

REGISTER_FUNCTION(atanh, eval_func(atanh_eval).
                         evalf_func(atanh_evalf).
                         derivative_func(atanh_deriv).
                         series_func(atanh_series).
                         conjugate_func(atanh_conjugate).
			 set_name("arctanh"));

//////////
// inverse hyperbolic cotangent (trigonometric function)
//////////

static ex acoth_evalf(const ex & x, PyObject* parent)
{
        if (is_exactly_a<numeric>(x))
                return acoth(ex_to<numeric>(x));

        return acoth(x).hold();
}

static ex acoth_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                // acoth(0) -> i*pi/2
                if (x.is_zero())
                        return _ex1_2*Pi*I;
                // acoth(1) -> oo
                if (x.is_equal(_ex1))
                        return Infinity;
                // acoth(-1) -> -oo
                if (x.is_equal(_ex_1))
                        return NegInfinity;
                //acoth(float) -> float 
                if (!x.info(info_flags::crational))
                        return acoth(ex_to<numeric>(x));
                // acoth() is odd
                if (x.info(info_flags::negative))
                        return -acoth(-x);
        }
       
        // acoth(oo) -> 0
        // acoth(-oo) -> 0
        // acoth(UnsignedInfinity) -> error
        if (x.info(info_flags::infinity)) {
                if (x.is_equal(Infinity) || x.is_equal(NegInfinity))
                        return _ex0;
                // x is UnsignedInfinity
                throw (std::runtime_error("arccoth_eval(): arccoth(unsigned_infinity) encountered"));
        }
        
        return acoth(x).hold();
}

static ex acoth_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // acoth(x) -> (1/2)*(ln(1 + 1/x) - ln(1 - 1/x))
        // d/dx acoth(x) -> 1/(1-x^2)
        return power(_ex1-power(x, _ex2), _ex_1);
}

static ex acoth_conjugate(const ex & x)
{
        // conjugate(acoth(x))==acoth(conjugate(x)) unless on the branch cuts which
        // run along the real axis inside the interval [-1, +1].
        if (is_exactly_a<numeric>(x) &&
            (!x.imag_part().is_zero() || (x < *_num_1_p && x > *_num1_p))) {
                return acoth(x.conjugate());
        }
        return conjugate_function(acoth(x)).hold();
}

REGISTER_FUNCTION(acoth, eval_func(acoth_eval).
                         evalf_func(acoth_evalf).
                         derivative_func(acoth_deriv).
                         conjugate_func(acoth_conjugate).
                         set_name("arccoth"));

} // namespace GiNaC
