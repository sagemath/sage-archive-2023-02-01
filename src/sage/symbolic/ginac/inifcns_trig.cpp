/** @file inifcns_trans.cpp
 *
 *  Implementation of trigonometric functions. */

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
#include "symbol.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "pseries.h"
#include "utils.h"


namespace GiNaC {


/* In Sage all the arc trig functions are printed with "arc" instead
   of "a" at the beginning.   This is for consistency with other
   computer algebra systems.   These print methods are registered
   below with each of the corresponding inverse trig function. */


//////////
// sine (trigonometric function)
//////////

static ex sin_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return sin(ex_to<numeric>(x));
	
	return sin(x).hold();
}

static ex sin_eval(const ex & x)
{
	// sin(oo) -> error
	// This should be before the tests below, since multiplying infinity
	// with other values raises runtime_errors
	if (x.info(info_flags::infinity)) {
		throw (std::runtime_error("sin_eval(): sin(infinity) encountered"));
	}

	// sin(n/d*Pi) -> { all known radicals with nesting depth 2 }
	const ex SixtyExOverPi = _ex60*x/Pi;
	ex sign = _ex1;
	if (is_exactly_a<numeric>(SixtyExOverPi)
                and SixtyExOverPi.info(info_flags::integer)) {
		numeric z = mod(ex_to<numeric>(SixtyExOverPi),*_num120_p);
		if (z>=*_num60_p) {
			// wrap to interval [0, Pi)
			z -= *_num60_p;
			sign = _ex_1;
		}
		if (z>*_num30_p) {
			// wrap to interval [0, Pi/2)
			z = *_num60_p-z;
		}
                // Not included were n*Pi/15 which has 3-level-deep roots
		if (z.is_equal(*_num0_p))  // sin(0)       -> 0
			return _ex0;
		if (z.is_equal(*_num2_p))  // sin(Pi/30)   -> -1/8*sqrt(5) + 1/2*sqrt(-3/8*sqrt(5) + 15/8) - 1/8
			return sign*(_ex_1*(sqrt(_ex5)+_ex1)/_ex8 +
                                _ex1_4*sqrt(_ex1_2*(_ex15-_ex3*sqrt(_ex5))));
		if (z.is_equal(*_num5_p))  // sin(Pi/12)   -> sqrt(6)/4*(1-sqrt(3)/3)
			return sign*_ex1_4*sqrt(_ex6)*(_ex1+_ex_1_3*sqrt(_ex3));
		if (z.is_equal(*_num6_p))  // sin(Pi/10)   -> sqrt(5)/4-1/4
			return sign*(_ex1_4*sqrt(_ex5)+_ex_1_4);
		if (z.is_equal(*_num10_p)) // sin(Pi/6)    -> 1/2
			return sign*_ex1_2;
		if (z.is_equal(*_num12_p)) // sin(Pi/5)    -> 1/4*sqrt(10-2*sqrt(5))
			return sign*_ex1_4*sqrt(_ex10-_ex2*sqrt(_ex5));
		if (z.is_equal(*_num14_p))  // sin(7*Pi/30)   -> -1/8*sqrt(5) + 1/2*sqrt(3/8*sqrt(5) + 15/8) + 1/8
			return sign*((_ex1-sqrt(_ex5))/_ex8 +
                                _ex1_4*sqrt(_ex1_2*(_ex15+_ex3*sqrt(_ex5))));
		if (z.is_equal(*_num15_p)) // sin(Pi/4)    -> sqrt(2)/2
			return sign*_ex1_2*sqrt(_ex2);
		if (z.is_equal(*_num18_p)) // sin(3/10*Pi) -> sqrt(5)/4+1/4
			return sign*(_ex1_4*sqrt(_ex5)+_ex1_4);
		if (z.is_equal(*_num20_p)) // sin(Pi/3)    -> sqrt(3)/2
			return sign*_ex1_2*sqrt(_ex3);
		if (z.is_equal(*_num22_p))  // sin(11*Pi/30)   -> 1/8*sqrt(5) + 1/2*sqrt(-3/8*sqrt(5) + 15/8) + 1/8
			return sign*((sqrt(_ex5)+_ex1)/_ex8 +
                                _ex1_4*sqrt(_ex1_2*(_ex15-_ex3*sqrt(_ex5))));
		if (z.is_equal(*_num24_p)) // sin(2*Pi/5)    -> 1/4*sqrt(10+2*sqrt(5))
			return sign*_ex1_4*sqrt(_ex10+_ex2*sqrt(_ex5));
		if (z.is_equal(*_num25_p)) // sin(5/12*Pi) -> sqrt(6)/4*(1+sqrt(3)/3)
			return sign*_ex1_4*sqrt(_ex6)*(_ex1+_ex1_3*sqrt(_ex3));
		if (z.is_equal(*_num26_p))  // sin(13*Pi/30)   -> 1/8*sqrt(5) + 1/2*sqrt(3/8*sqrt(5) + 15/8) - 1/8
			return sign*((sqrt(_ex5)-_ex1)/_ex8 +
                                _ex1_4*sqrt(_ex1_2*(_ex15+_ex3*sqrt(_ex5))));
		if (z.is_equal(*_num30_p)) // sin(Pi/2)    -> 1
			return sign;
	}

	const ex TwentyforExOverPi = _ex24*x/Pi;
        sign = _ex1;
	if (is_exactly_a<numeric>(TwentyforExOverPi)
                and TwentyforExOverPi.info(info_flags::integer)) {
		numeric z = mod(ex_to<numeric>(TwentyforExOverPi),*_num48_p);
		if (z>=*_num24_p) {
			// wrap to interval [0, Pi)
			z -= *_num24_p;
			sign = _ex_1;
		}
		if (z>*_num12_p) {
			// wrap to interval [0, Pi/2)
			z = *_num24_p-z;
		}
		if (z.is_equal(*_num1_p))  // sin(Pi/24) -> 1/4*sqrt(8-2*sqrt(6)-2*sqrt(2))
			return sign*(_ex1_4*sqrt(_ex8 - _ex2*sqrt(_ex6) - _ex2*sqrt(_ex2)));
		if (z.is_equal(*_num3_p))  // sin(Pi/8) -> 1/2*sqrt(-sqrt(2) + 2)
			return sign*(_ex1_2*sqrt(_ex2-sqrt(_ex2)));
		if (z.is_equal(*_num5_p))  // sin(5*Pi/24) -> 1/4*sqrt(8-2*sqrt(6)+2*sqrt(2))
			return sign*(_ex1_4*sqrt(_ex8 - _ex2*sqrt(_ex6) + _ex2*sqrt(_ex2)));
		if (z.is_equal(*_num7_p))  // sin(7*Pi/24) -> 1/4*sqrt(8+2*sqrt(6)-2*sqrt(2))
			return sign *(_ex1_4*sqrt(_ex8 + _ex2*sqrt(_ex6) - _ex2*sqrt(_ex2)));
		if (z.is_equal(*_num9_p))  // sin(3*Pi/8) -> 1/2*sqrt(-sqrt(2) + 2)
			return sign*(_ex1_2*sqrt(_ex2+sqrt(_ex2)));
		if (z.is_equal(*_num11_p))  // sin(11*Pi/24) -> 1/4*sqrt(8+2*sqrt(6)+2*sqrt(2))
			return sign *(_ex1_4*sqrt(_ex8 + _ex2*sqrt(_ex6) + _ex2*sqrt(_ex2)));
	}

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// sin(asin(x)) -> x
		if (is_ex_the_function(x, asin))
			return t;

		// sin(acos(x)) -> sqrt(1-x^2)
		if (is_ex_the_function(x, acos))
			return sqrt(_ex1-power(t,_ex2));

		// sin(atan(x)) -> x/sqrt(1+x^2)
		if (is_ex_the_function(x, atan))
			return t*power(_ex1+power(t,_ex2),_ex_1_2);
	}

	// try something
	const ex ExOverPi = x/Pi;
	if (ExOverPi.info(info_flags::integer))
		return _ex0;
	
	// sin(float) -> float
        if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
		return sin(ex_to<numeric>(x));

	// sin() is odd
	if (x.info(info_flags::negative))
		return -sin(-x);
	
	return sin(x).hold();
}

static ex sin_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx sin(x) -> cos(x)
	return cos(x);
}

static ex sin_real_part(const ex & x)
{
	return cosh(GiNaC::imag_part(x))*sin(GiNaC::real_part(x));
}

static ex sin_imag_part(const ex & x)
{
	return sinh(GiNaC::imag_part(x))*cos(GiNaC::real_part(x));
}

static ex sin_conjugate(const ex & x)
{
	// conjugate(sin(x))==sin(conjugate(x))
	return sin(x.conjugate());
}

REGISTER_FUNCTION(sin, eval_func(sin_eval).
                       evalf_func(sin_evalf).
                       derivative_func(sin_deriv).
                       real_part_func(sin_real_part).
                       imag_part_func(sin_imag_part).
                       conjugate_func(sin_conjugate).
                       latex_name("\\sin"));

//////////
// cosine (trigonometric function)
//////////

static ex cos_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return cos(ex_to<numeric>(x));
	
	return cos(x).hold();
}

static ex cos_eval(const ex & x)
{
	// cos(oo) -> error
	// This should be before the tests below, since multiplying infinity
	// with other values raises runtime_errors
	if (x.info(info_flags::infinity)) {
		throw (std::runtime_error("cos_eval(): cos(infinity) encountered"));
	}

	// cos(n/d*Pi) -> { all known radicals with nesting depth 2 }
	const ex SixtyExOverPi = _ex60*x/Pi;
	ex sign = _ex1;
	if (is_exactly_a<numeric>(SixtyExOverPi)
                and SixtyExOverPi.info(info_flags::integer)) {
		numeric z = mod(ex_to<numeric>(SixtyExOverPi),*_num120_p);
		if (z>=*_num60_p) {
			// wrap to interval [0, Pi)
			z = *_num120_p-z;
		}
		if (z>=*_num30_p) {
			// wrap to interval [0, Pi/2)
			z = *_num60_p-z;
			sign = _ex_1;
		}
		if (z.is_equal(*_num0_p))  // cos(0)       -> 1
			return sign;
		if (z.is_equal(*_num4_p))  // cos(Pi/15)   -> 1/8*sqrt(5) + 1/2*sqrt(3/8*sqrt(5) + 15/8) - 1/8
			return sign*((sqrt(_ex5)-_ex1)/_ex8 + _ex1_4*sqrt((_ex15+_ex3*sqrt(_ex5))/_ex2));
		if (z.is_equal(*_num5_p))  // cos(Pi/12)   -> sqrt(6)/4*(1+sqrt(3)/3)
			return sign*_ex1_4*sqrt(_ex6)*(_ex1+_ex1_3*sqrt(_ex3));
		if (z.is_equal(*_num6_p))  // cos(Pi/10)   -> 1/2*sqrt(1/2*sqrt(5) + 5/2)
			return sign*_ex1_2*sqrt(_ex1_2*(_ex5+sqrt(_ex5)));
		if (z.is_equal(*_num8_p))  // cos(2*Pi/15)   -> 1/8*sqrt(5) + 1/2*sqrt(-3/8*sqrt(5) + 15/8) + 1/8
			return sign*((sqrt(_ex5)+_ex1)/_ex8 + _ex1_4*sqrt((_ex15-_ex3*sqrt(_ex5))/_ex2));
		if (z.is_equal(*_num10_p)) // cos(Pi/6)    -> sqrt(3)/2
			return sign*_ex1_2*sqrt(_ex3);
		if (z.is_equal(*_num12_p)) // cos(Pi/5)    -> sqrt(5)/4+1/4
			return sign*(_ex1_4*sqrt(_ex5)+_ex1_4);
		if (z.is_equal(*_num15_p)) // cos(Pi/4)    -> sqrt(2)/2
			return sign*_ex1_2*sqrt(_ex2);
		if (z.is_equal(*_num16_p))  // cos(4*Pi/15)   -> -1/8*sqrt(5) + 1/2*sqrt(3/8*sqrt(5) + 15/8) + 1/8
			return sign*((_ex1-sqrt(_ex5))/_ex8 + _ex1_4*sqrt((_ex15+_ex3*sqrt(_ex5))/_ex2));
		if (z.is_equal(*_num18_p))  // cos(3*Pi/10)   -> 1/2*sqrt(-1/2*sqrt(5) + 5/2)
			return sign*_ex1_2*sqrt(_ex1_2*(_ex5-sqrt(_ex5)));
		if (z.is_equal(*_num20_p)) // cos(Pi/3)    -> 1/2
			return sign*_ex1_2;
		if (z.is_equal(*_num24_p)) // cos(2/5*Pi)  -> sqrt(5)/4-1/4x
			return sign*(_ex1_4*sqrt(_ex5)+_ex_1_4);
		if (z.is_equal(*_num25_p)) // cos(5/12*Pi) -> sqrt(6)/4*(1-sqrt(3)/3)
			return sign*_ex1_4*sqrt(_ex6)*(_ex1+_ex_1_3*sqrt(_ex3));
		if (z.is_equal(*_num28_p))  // cos(7*Pi/15)   -> -1/8*sqrt(5) + 1/2*sqrt(-3/8*sqrt(5) + 15/8) - 1/8
			return sign*(_ex_1*(_ex1+sqrt(_ex5))/_ex8 + _ex1_4*sqrt((_ex15-_ex3*sqrt(_ex5))/_ex2));
		if (z.is_equal(*_num30_p)) // cos(Pi/2)    -> 0
			return _ex0;
	}

	const ex TwentyforExOverPi = _ex24*x/Pi;
        sign = _ex1;
	if (is_exactly_a<numeric>(TwentyforExOverPi)
                        and TwentyforExOverPi.info(info_flags::integer)) {
		numeric z = mod(ex_to<numeric>(TwentyforExOverPi),*_num48_p);
		if (z>=*_num24_p) {
			// wrap to interval [0, Pi)
			z -= *_num48_p;
		}
		if (z>*_num12_p) {
			// wrap to interval [0, Pi/2)
			z = *_num24_p-z;
			sign = _ex_1;
		}
		if (z.is_equal(*_num1_p))  // cos(Pi/24) -> 1/4*sqrt(8+2*sqrt(6)+2*sqrt(2))
			return sign*(_ex1_4*sqrt(_ex8 + _ex2*sqrt(_ex6) + _ex2*sqrt(_ex2)));
		if (z.is_equal(*_num3_p))  // cos(Pi/8) -> 1/2*sqrt(sqrt(2) + 2)
			return sign*(_ex1_2*sqrt(_ex2+sqrt(_ex2)));
		if (z.is_equal(*_num5_p))  // cos(5*Pi/24) -> 1/4*sqrt(8+2*sqrt(6)-2*sqrt(2))
			return sign*(_ex1_4*sqrt(_ex8 + _ex2*sqrt(_ex6) - _ex2*sqrt(_ex2)));
		if (z.is_equal(*_num7_p))  // cos(7*Pi/24) -> 1/4*sqrt(8-2*sqrt(6)+2*sqrt(2))
			return sign*(_ex1_4*sqrt(_ex8 - _ex2*sqrt(_ex6) + _ex2*sqrt(_ex2)));
		if (z.is_equal(*_num9_p))  // cos(3*Pi/8) -> 1/2*sqrt(-sqrt(2) + 2)
			return sign*(_ex1_2*sqrt(_ex2-sqrt(_ex2)));
		if (z.is_equal(*_num11_p))  // cos(11*Pi/24) -> 1/4*sqrt(8-2*sqrt(6)-2*sqrt(2))
			return sign*(_ex1_4*sqrt(_ex8 - _ex2*sqrt(_ex6) - _ex2*sqrt(_ex2)));
	}

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// cos(acos(x)) -> x
		if (is_ex_the_function(x, acos))
			return t;

		// cos(asin(x)) -> sqrt(1-x^2)
		if (is_ex_the_function(x, asin))
			return sqrt(_ex1-power(t,_ex2));

		// cos(atan(x)) -> 1/sqrt(1+x^2)
		if (is_ex_the_function(x, atan))
			return power(_ex1+power(t,_ex2),_ex_1_2);
	}

	// cos(integer*pi) --> (-1)^integer
	const ex ExOverPi = x/Pi;
	if (ExOverPi.info(info_flags::integer))
		return pow(_ex_1, ExOverPi);
	
	// cos(float) -> float
	if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
		return cos(ex_to<numeric>(x));
	
	// cos() is even
	if (x.info(info_flags::negative))
		return cos(-x);
	
	return cos(x).hold();
}

static ex cos_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx cos(x) -> -sin(x)
	return -sin(x);
}

static ex cos_real_part(const ex & x)
{
	return cosh(GiNaC::imag_part(x))*cos(GiNaC::real_part(x));
}

static ex cos_imag_part(const ex & x)
{
	return -sinh(GiNaC::imag_part(x))*sin(GiNaC::real_part(x));
}

static ex cos_conjugate(const ex & x)
{
	// conjugate(cos(x))==cos(conjugate(x))
	return cos(x.conjugate());
}

REGISTER_FUNCTION(cos, eval_func(cos_eval).
                       evalf_func(cos_evalf).
                       derivative_func(cos_deriv).
                       real_part_func(cos_real_part).
                       imag_part_func(cos_imag_part).
                       conjugate_func(cos_conjugate).
                       latex_name("\\cos"));

//////////
// tangent (trigonometric function)
//////////

static ex tan_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return tan(ex_to<numeric>(x));
	
	return tan(x).hold();
}

static ex tan_eval(const ex & x)
{
	// tan(oo) -> error
	// This should be before the tests below, since multiplying infinity
	// with other values raises runtime_errors
	if (x.info(info_flags::infinity)) {
		throw (std::runtime_error("tan_eval(): tan(infinity) encountered"));
	}

	// tan(n/d*Pi) -> { all known non-nested radicals }
	const ex SixtyExOverPi = _ex60*x/Pi;
	ex sign = _ex1;
	if (is_exactly_a<numeric>(SixtyExOverPi)
                and SixtyExOverPi.info(info_flags::integer)) {
		numeric z = mod(ex_to<numeric>(SixtyExOverPi),*_num60_p);
		if (z>=*_num60_p) {
			// wrap to interval [0, Pi)
			z -= *_num60_p;
		}
		if (z>=*_num30_p) {
			// wrap to interval [0, Pi/2)
			z = *_num60_p-z;
			sign = _ex_1;
		}
		if (z.is_equal(*_num0_p))  // tan(0)       -> 0
			return _ex0;
		if (z.is_equal(*_num3_p)) // tan(Pi/20)    -> sqrt(5) - 1/2*sqrt(8*sqrt(5) + 20) + 1
			return sign*(sqrt(_ex5)+_ex1-sqrt(_ex20+_ex8*sqrt(_ex5))/_ex2);
		if (z.is_equal(*_num5_p))  // tan(Pi/12)   -> 2-sqrt(3)
			return sign*(_ex2-sqrt(_ex3));
		if (z.is_equal(*_num6_p)) // tan(Pi/10)    -> sqrt(1-2/5*sqrt(5))
			return sign*sqrt(_ex1-_ex2/_ex5*sqrt(_ex5));
		if (z.is_equal(*_num9_p)) // tan(3*Pi/20)    -> sqrt(5) - 1/2*sqrt(-8*sqrt(5) + 20) - 1
			return sign*(sqrt(_ex5)-_ex1-sqrt(_ex20-_ex8*sqrt(_ex5))/_ex2);
		if (z.is_equal(*_num10_p)) // tan(Pi/6)    -> sqrt(3)/3
			return sign*_ex1_3*sqrt(_ex3);
		if (z.is_equal(*_num12_p)) // tan(Pi/5)    -> sqrt(5-2*sqrt(5))
			return sign*sqrt(_ex5-_ex2*sqrt(_ex5));
		if (z.is_equal(*_num15_p)) // tan(Pi/4)    -> 1
			return sign;
		if (z.is_equal(*_num18_p)) // tan(3*Pi/10)    -> sqrt(1+2/5*sqrt(5))
			return sign*sqrt(_ex1+_ex2/_ex5*sqrt(_ex5));
		if (z.is_equal(*_num20_p)) // tan(Pi/3)    -> sqrt(3)
			return sign*sqrt(_ex3);
		if (z.is_equal(*_num21_p)) // tan(7*Pi/20)    -> sqrt(5) + 1/2*sqrt(-8*sqrt(5) + 20) - 1
			return sign*(sqrt(_ex5)-_ex1+sqrt(_ex20-_ex8*sqrt(_ex5))/_ex2);
		if (z.is_equal(*_num24_p)) // tan(2*Pi/5)    -> sqrt(5+2*sqrt(5))
			return sign*sqrt(_ex5+_ex2*sqrt(_ex5));
		if (z.is_equal(*_num25_p)) // tan(5/12*Pi) -> 2+sqrt(3)
			return sign*(sqrt(_ex3)+_ex2);
		if (z.is_equal(*_num27_p)) // tan(9*Pi/20)    -> sqrt(5) + 1/2*sqrt(8*sqrt(5) + 20) + 1
			return sign*(sqrt(_ex5)+_ex1+sqrt(_ex20+_ex8*sqrt(_ex5))/_ex2);
		if (z.is_equal(*_num30_p)) // tan(Pi/2)    -> infinity
			//throw (pole_error("tan_eval(): simple pole",1));
			return UnsignedInfinity;
	}

	const ex FortyeightExOverPi = _ex48*x/Pi;
        sign = _ex1;
	if (is_exactly_a<numeric>(FortyeightExOverPi)
                        and FortyeightExOverPi.info(info_flags::integer)) {
		numeric z = mod(ex_to<numeric>(FortyeightExOverPi),*_num48_p);
		if (z>=*_num48_p) {
			// wrap to interval [0, Pi)
			z -= *_num48_p;
		}
		if (z>*_num24_p) {
			// wrap to interval [0, Pi/2)
			z = *_num48_p-z;
			sign = _ex_1;
		}
		if (z.is_equal(*_num2_p))  // tan(Pi/24) -> sqrt(6) - sqrt(-2*sqrt(6) + 5) - 2
			return sign *(_ex_2+sqrt(_ex6)+sqrt(_ex2)-sqrt(_ex3));
		if (z.is_equal(*_num3_p))  // tan(Pi/16) -> -sqrt(2) + sqrt(2*sqrt(2) + 4) - 1
			return sign*(_ex_1-sqrt(_ex2)+sqrt(_ex2*sqrt(_ex2)+_ex4));
		if (z.is_equal(*_num6_p))  // tan(Pi/8) -> sqrt(2)-1
			return sign*(sqrt(_ex2)-_ex1);
		if (z.is_equal(*_num9_p))  // tan(3*Pi/16) -> -sqrt(2) + sqrt(-2*sqrt(2) + 4) + 1
			return sign*(_ex1-sqrt(_ex2)+sqrt(_ex4-_ex2*sqrt(_ex2)));
		if (z.is_equal(*_num10_p))  // tan(5*Pi/24) -> sqrt(6) + sqrt(-2*sqrt(6) + 5) - 2
			return sign*(_ex_2+sqrt(_ex6)-sqrt(_ex2)+sqrt(_ex3));
		if (z.is_equal(*_num14_p))  // tan(7*Pi/24) -> sqrt(6) - sqrt(2*sqrt(6) + 5) + 2
			return sign*(_ex2+sqrt(_ex6)-sqrt(_ex2)-sqrt(_ex3));
		if (z.is_equal(*_num15_p))  // tan(5*Pi/16) -> sqrt(2) + sqrt(-2*sqrt(2) + 4) + 1
			return sign*(_ex_1+sqrt(_ex2)+sqrt(_ex4-_ex2*sqrt(_ex2)));
		if (z.is_equal(*_num18_p))  // tan(3*Pi/8) -> sqrt(2)+1
			return sign*(sqrt(_ex2)+_ex1);
		if (z.is_equal(*_num21_p))  // tan(7*Pi/16) -> sqrt(2) + sqrt(2*sqrt(2) + 4) + 1
			return sign*(_ex1+sqrt(_ex2)+sqrt(_ex2*sqrt(_ex2)+_ex4));
		if (z.is_equal(*_num22_p))  // tan(11*Pi/24) -> sqrt(6) + sqrt(2*sqrt(6) + 5) + 2
			return sign*(_ex2+sqrt(_ex6)+sqrt(_ex2)+sqrt(_ex3));
	}

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// tan(atan(x)) -> x
		if (is_ex_the_function(x, atan))
			return t;

		// tan(asin(x)) -> x/sqrt(1+x^2)
		if (is_ex_the_function(x, asin))
			return t*power(_ex1-power(t,_ex2),_ex_1_2);

		// tan(acos(x)) -> sqrt(1-x^2)/x
		if (is_ex_the_function(x, acos))
			return power(t,_ex_1)*sqrt(_ex1-power(t,_ex2));
	}
	
	// tan(integer*pi) -> 0
	const ex ExOverPi = x/Pi;
	if (ExOverPi.info(info_flags::integer))
		return _ex0;

	// tan(float) -> float
	if (x.info(info_flags::numeric) && !x.info(info_flags::crational)) {
		return tan(ex_to<numeric>(x));
	}
	
	// tan() is odd
	if (x.info(info_flags::negative))
		return -tan(-x);
	
	return tan(x).hold();
}

static ex tan_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx tan(x) -> 1+tan(x)^2;
	return (_ex1+power(tan(x),_ex2));
}

static ex tan_real_part(const ex & x)
{
	ex a = GiNaC::real_part(mul(x, _ex2));
	ex b = GiNaC::imag_part(mul(x, _ex2));
	return sin(a)/(cos(a) + cosh(b));
}

static ex tan_imag_part(const ex & x)
{
	ex a = GiNaC::real_part(mul(x, _ex2));
	ex b = GiNaC::imag_part(mul(x, _ex2));
	return sinh(b)/(cos(a) + cosh(b));
}

static ex tan_series(const ex &x,
                     const relational &rel,
                     int order,
                     unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole falls back to tan_deriv.
	// On a pole simply expand sin(x)/cos(x).
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (!(2*x_pt/Pi).info(info_flags::odd))
		throw do_taylor();  // caught by function::series()
	// if we got here we have to care for a simple pole
	return (sin(x)/cos(x)).series(rel, order, options);
}

static ex tan_conjugate(const ex & x)
{
	// conjugate(tan(x))==tan(conjugate(x))
	return tan(x.conjugate());
}

REGISTER_FUNCTION(tan, eval_func(tan_eval).
                       evalf_func(tan_evalf).
                       derivative_func(tan_deriv).
                       series_func(tan_series).
                       real_part_func(tan_real_part).
                       imag_part_func(tan_imag_part).
                       conjugate_func(tan_conjugate).
                       latex_name("\\tan"));

//////////
// cotangent (trigonometric function)
//////////

static ex cot_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return tan(ex_to<numeric>(x)).inverse();

	return cot(x).hold();
}

static ex cot_eval(const ex & x)
{
	// cot(0) -> error
	// This should be before the tests below, since multiplying infinity
	// with other values raises runtime_errors
	if (x.is_zero()) {
		throw (std::runtime_error("cotan_eval(): cot(0) encountered"));
	}

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// tan(atan(x)) -> x
		if (is_ex_the_function(x, acot))
			return t;
	}

	// cot(float) -> float
	if (x.info(info_flags::numeric) && !x.info(info_flags::crational)) {
		return tan(ex_to<numeric>(x)).inverse();
	}

        ex res = tan_eval(x);
	if (not is_ex_the_function(res, tan))
                return power(res, _ex_1);

	return cot(x).hold();
}

static ex cot_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx cot(x) -> -1-cot(x)^2;
	return (_ex_1-power(cot(x),_ex2));
}

// Ref.: http://dlmf.nist.gov/4.21.E40
static ex cot_real_part(const ex & x)
{
	ex a = GiNaC::real_part(mul(x, _ex2));
	ex b = GiNaC::imag_part(mul(x, _ex2));
	return sin(a)/(cosh(b) - cos(a));
}

static ex cot_imag_part(const ex & x)
{
	ex a = GiNaC::real_part(mul(x, _ex2));
	ex b = GiNaC::imag_part(mul(x, _ex2));
	return sinh(b)/(cosh(b) - cos(a));
}

static ex cot_series(const ex &x,
                     const relational &rel,
                     int order,
                     unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole falls back to tan_deriv.
	// On a pole simply expand cos(x)/sin(x).
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (!(2*x_pt/Pi).info(info_flags::even))
		throw do_taylor();  // caught by function::series()
	// if we got here we have to care for a simple pole
	return (cos(x)/sin(x)).series(rel, order, options);
}

static ex cot_conjugate(const ex & x)
{
	// conjugate(tan(x))==1/tan(conjugate(x))
	return power(tan(x.conjugate()), _ex_1);
}

REGISTER_FUNCTION(cot, eval_func(cot_eval).
                       evalf_func(cot_evalf).
                       derivative_func(cot_deriv).
                       series_func(cot_series).
                       real_part_func(cot_real_part).
                       imag_part_func(cot_imag_part).
                       conjugate_func(cot_conjugate).
                       latex_name("\\cot"));

//////////
// inverse sine (arc sine)
//////////

static ex asin_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return asin(ex_to<numeric>(x));
	
	return asin(x).hold();
}

static ex asin_eval(const ex & x)
{
	if (x.info(info_flags::numeric)) {

		// asin(0) -> 0
		if (x.is_zero())
			return x;

		// asin(1/2) -> Pi/6
		if (x.is_equal(_ex1_2))
			return numeric(1,6)*Pi;

		// asin(1) -> Pi/2
		if (x.is_equal(_ex1))
			return _ex1_2*Pi;

		// asin(-1/2) -> -Pi/6
		if (x.is_equal(_ex_1_2))
			return numeric(-1,6)*Pi;

		// asin(-1) -> -Pi/2
		if (x.is_equal(_ex_1))
			return _ex_1_2*Pi;

		// asin(float) -> float
		if (!x.info(info_flags::crational))
			return asin(ex_to<numeric>(x));

		// asin() is odd
		if (x.info(info_flags::negative))
			return -asin(-x);
	}
	
	// asin(oo) -> error
	// asin(UnsignedInfinity) -> UnsignedInfinity
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(UnsignedInfinity))
			return UnsignedInfinity;
		throw (std::runtime_error("arcsin_eval(): arcsin(infinity) encountered"));
	}

	return asin(x).hold();
}

static ex asin_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx asin(x) -> 1/sqrt(1-x^2)
	return power(1-power(x,_ex2),_ex_1_2);
}

static ex asin_conjugate(const ex & x)
{
	// conjugate(asin(x))==asin(conjugate(x)) unless on the branch cuts which
	// run along the real axis outside the interval [-1, +1].
	if (is_exactly_a<numeric>(x) &&
	    (!x.imag_part().is_zero() || (x > *_num_1_p && x < *_num1_p))) {
		return asin(x.conjugate());
	}
	return conjugate_function(asin(x)).hold();
}

REGISTER_FUNCTION(asin, eval_func(asin_eval).
                        evalf_func(asin_evalf).
                        derivative_func(asin_deriv).
                        conjugate_func(asin_conjugate).
			set_name("arcsin", "\\arcsin"));

//////////
// inverse cosine (arc cosine)
//////////

static ex acos_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return acos(ex_to<numeric>(x));
	
	return acos(x).hold();
}

static ex acos_eval(const ex & x)
{
	if (x.info(info_flags::numeric)) {

		// acos(1) -> 0
		if (x.is_equal(_ex1))
			return _ex0;

		// acos(1/2) -> Pi/3
		if (x.is_equal(_ex1_2))
			return _ex1_3*Pi;

		// acos(0) -> Pi/2
		if (x.is_zero())
			return _ex1_2*Pi;

		// acos(-1/2) -> 2/3*Pi
		if (x.is_equal(_ex_1_2))
			return numeric(2,3)*Pi;

		// acos(-1) -> Pi
		if (x.is_equal(_ex_1))
			return Pi;

		// acos(float) -> float
		if (!x.info(info_flags::crational))
			return acos(ex_to<numeric>(x));

		// acos(-x) -> Pi-acos(x)
		if (x.info(info_flags::negative))
			return Pi-acos(-x);
	}
	
	// acos(oo) -> error
	// acos(UnsignedInfinity) -> UnsignedInfinity
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(UnsignedInfinity))
			return UnsignedInfinity;
		throw (std::runtime_error("arccos_eval(): arccos(infinity) encountered"));
	}
	return acos(x).hold();
}

static ex acos_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx acos(x) -> -1/sqrt(1-x^2)
	return -power(1-power(x,_ex2),_ex_1_2);
}

static ex acos_conjugate(const ex & x)
{
	// conjugate(acos(x))==acos(conjugate(x)) unless on the branch cuts which
	// run along the real axis outside the interval [-1, +1].
	if (is_exactly_a<numeric>(x) &&
	    (!x.imag_part().is_zero() || (x > *_num_1_p && x < *_num1_p))) {
		return acos(x.conjugate());
	}
	return conjugate_function(acos(x)).hold();
}

 
REGISTER_FUNCTION(acos, eval_func(acos_eval).
                        evalf_func(acos_evalf).
                        derivative_func(acos_deriv).
                        conjugate_func(acos_conjugate).
			set_name("arccos", "\\arccos"));

//////////
// inverse tangent (arc tangent)
//////////

static ex atan_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return atan(ex_to<numeric>(x));
	
	return atan(x).hold();
}

static ex atan_eval(const ex & x)
{
	if (x.info(info_flags::numeric)) {

		// atan(0) -> 0
		if (x.is_zero())
			return _ex0;

		// atan(1) -> Pi/4
		if (x.is_equal(_ex1))
			return _ex1_4*Pi;

		// atan(-1) -> -Pi/4
		if (x.is_equal(_ex_1))
			return _ex_1_4*Pi;

		if (x.is_equal(I) || x.is_equal(-I))
			throw (pole_error("atan_eval(): logarithmic pole",0));

		// atan(float) -> float
		if (!x.info(info_flags::crational))
			return atan(ex_to<numeric>(x));

		// atan() is odd
		if (x.info(info_flags::negative))
			return -atan(-x);
	}
	
	// arctan(oo) -> Pi/2
	// arctan(-oo) -> -Pi/2
	// arctan(UnsignedInfinity) -> error
	if (is_exactly_a<infinity>(x)) {
	        infinity xinf = ex_to<infinity>(x);
		if (xinf.is_plus_infinity())
		        return _ex1_2*Pi;
		if (xinf.is_minus_infinity())
		        return _ex_1_2*Pi;
		// x is UnsignedInfinity
		throw (std::runtime_error("arctan_eval(): arctan(unsigned_infinity) encountered"));
	}
		
	return atan(x).hold();
}

static ex atan_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx atan(x) -> 1/(1+x^2)
	return power(_ex1+power(x,_ex2), _ex_1);
}

static ex atan_series(const ex &arg,
                      const relational &rel,
                      int order,
                      unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole or cut falls back to atan_deriv.
	// There are two branch cuts, one runnig from I up the imaginary axis and
	// one running from -I down the imaginary axis.  The points I and -I are
	// poles.
	// On the branch cuts and the poles series expand
	//     (log(1+I*x)-log(1-I*x))/(2*I)
	// instead.
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (!(I*arg_pt).info(info_flags::real))
		throw do_taylor();     // Re(x) != 0
	if ((I*arg_pt).info(info_flags::real) && abs(I*arg_pt)<_ex1)
		throw do_taylor();     // Re(x) == 0, but abs(x)<1
	// care for the poles, using the defining formula for atan()...
	if (arg_pt.is_equal(I) || arg_pt.is_equal(-I))
		return ((log(1+I*arg)-log(1-I*arg))/(2*I)).series(rel, order, options);
	if (!(options & series_options::suppress_branchcut)) {
		// method:
		// This is the branch cut: assemble the primitive series manually and
		// then add the corresponding complex step function.
		const symbol &s = ex_to<symbol>(rel.lhs());
		const ex &point = rel.rhs();
		const symbol foo;
		const ex replarg = series(atan(arg), s==foo, order).subs(foo==point, subs_options::no_pattern);
		ex Order0correction = replarg.op(0)+csgn(arg)*Pi*_ex_1_2;
		if ((I*arg_pt)<_ex0)
			Order0correction += log((I*arg_pt+_ex_1)/(I*arg_pt+_ex1))*I*_ex_1_2;
		else
			Order0correction += log((I*arg_pt+_ex1)/(I*arg_pt+_ex_1))*I*_ex1_2;
		epvector seq;
		seq.push_back(expair(Order0correction, _ex0));
		seq.push_back(expair(Order(_ex1), order));
		return series(replarg - pseries(rel, seq), rel, order);
	}
	throw do_taylor();
}

static ex atan_conjugate(const ex & x)
{
	// conjugate(atan(x))==atan(conjugate(x)) unless on the branch cuts which
	// run along the imaginary axis outside the interval [-I, +I].
	if (x.info(info_flags::real))
		return atan(x);
	if (is_exactly_a<numeric>(x)) {
		const numeric x_re = ex_to<numeric>(x.real_part());
		const numeric x_im = ex_to<numeric>(x.imag_part());
		if (!x_re.is_zero() ||
		    (x_im > *_num_1_p && x_im < *_num1_p))
			return atan(x.conjugate());
	}
	return conjugate_function(atan(x)).hold();
}

REGISTER_FUNCTION(atan, eval_func(atan_eval).
                        evalf_func(atan_evalf).
                        derivative_func(atan_deriv).
                        series_func(atan_series).
                        conjugate_func(atan_conjugate).
			set_name("arctan", "\\arctan"));

//////////
// inverse tangent (atan2(y,x))
//////////

static ex atan2_evalf(const ex &y, const ex &x, PyObject* parent)
{
	if (is_exactly_a<numeric>(y) && is_exactly_a<numeric>(x))
		return atan(ex_to<numeric>(y), ex_to<numeric>(x));
	
	return atan2(y, x).hold();
}

static ex atan2_eval(const ex & y, const ex & x)
{
	if (y.is_zero()) {

		// atan2(0, 0) -> undefined
		if (x.is_zero())
			throw (std::runtime_error("arctan2_eval(): arctan2(0,0) encountered"));

		// atan2(0, x), x real and positive -> 0
		if (x.info(info_flags::positive))
			return _ex0;

		// atan2(0, x), x real and negative -> Pi
		if (x.info(info_flags::negative))
			return Pi;
	}

	if (x.is_zero()) {

		// atan2(y, 0), y real and positive -> Pi/2
		if (y.info(info_flags::positive))
			return _ex1_2*Pi;

		// atan2(y, 0), y real and negative -> -Pi/2
		if (y.info(info_flags::negative))
			return _ex_1_2*Pi;
	}

	if (y.is_equal(x)) {

		// atan2(y, y), y real and positive -> Pi/4
		if (y.info(info_flags::positive))
			return _ex1_4*Pi;

		// atan2(y, y), y real and negative -> -3/4*Pi
		if (y.info(info_flags::negative))
			return numeric(-3, 4)*Pi;
	}

	if (y.is_equal(-x)) {

		// atan2(y, -y), y real and positive -> 3*Pi/4
		if (y.info(info_flags::positive))
			return numeric(3, 4)*Pi;

		// atan2(y, -y), y real and negative -> -Pi/4
		if (y.info(info_flags::negative))
			return _ex_1_4*Pi;
	}

	// atan2(float, float) -> float
	if (is_a<numeric>(y) && !y.info(info_flags::crational) &&
	    is_a<numeric>(x) && !x.info(info_flags::crational))
		return atan(ex_to<numeric>(y), ex_to<numeric>(x));

	// handle infinities
	if (is_a<infinity>(x) || is_a<infinity>(y)) {
		if (is_a<infinity>(x) && ex_to<infinity>(x).is_unsigned_infinity())
			throw (std::runtime_error("arctan2_eval(): arctan2(unsigned_infinity, x) encountered"));
		if (is_a<infinity>(y) && ex_to<infinity>(y).is_unsigned_infinity())
			throw (std::runtime_error("arctan2_eval(): arctan2(x, unsigned_infinity) encountered"));

		if (is_a<infinity>(x) && is_a<infinity>(y)) 
			return atan2_eval(ex_to<infinity>(x).get_direction(), 
					  ex_to<infinity>(y).get_direction());

		if (is_a<infinity>(x)) 
			return atan2_eval(ex_to<infinity>(x).get_direction(), 0);
		if (is_a<infinity>(y)) 
			return atan2_eval(0, ex_to<infinity>(y).get_direction());
	}

	// atan2(real, real) -> atan(y/x) +/- Pi
	if (y.info(info_flags::real) && x.info(info_flags::real)) {
		if (x.info(info_flags::positive))
			return atan(y/x);

		if (x.info(info_flags::negative)) {
			if (y.info(info_flags::positive))
				return atan(y/x)+Pi;
			if (y.info(info_flags::negative))
				return atan(y/x)-Pi;
		}
	}
		
	return atan2(y, x).hold();
}    

static ex atan2_deriv(const ex & y, const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param<2);
	
	if (deriv_param==0) {
		// d/dy atan2(y,x)
		return x*power(power(x,_ex2)+power(y,_ex2),_ex_1);
	}
	// d/dx atan2(y,x)
	return -y*power(power(x,_ex2)+power(y,_ex2),_ex_1);
}

REGISTER_FUNCTION(atan2, eval_func(atan2_eval).
                         evalf_func(atan2_evalf).
                         derivative_func(atan2_deriv).
			 set_name("arctan2"));

//////////
// inverse cotangent (arc cotangent)
//////////

static ex acot_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return atan(ex_to<numeric>(x).inverse());

	return acot(x).hold();
}

static ex acot_eval(const ex & x)
{
	if (x.info(info_flags::numeric)) {

		if (x.is_zero())
			return _ex1_2*Pi;

		if (x.is_equal(_ex1))
			return _ex1_4*Pi;

		if (x.is_equal(_ex_1))
			return _ex_1_4*Pi;

		if (x.is_equal(I) || x.is_equal(-I))
			throw (pole_error("acot_eval(): logarithmic pole",0));

                if (!x.info(info_flags::crational))
                        return atan(ex_to<numeric>(x).inverse());

		if (x.info(info_flags::negative))
			return -acot(-x);
	}

	if (x.info(info_flags::infinity)) {
		return _ex0;
	}

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);
		if (is_ex_the_function(x, cot))
			return t;
	}
	return acot(x).hold();
}

static ex acot_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	return mul(power(_ex1+power(x,_ex2), _ex_1), _ex_1);
}

static ex acot_series(const ex &arg,
                      const relational &rel,
                      int order,
                      unsigned options)
{
        return _ex1_2*Pi - atan_series(arg, rel, order, options);
}

REGISTER_FUNCTION(acot, eval_func(acot_eval).
                        evalf_func(acot_evalf).
                        derivative_func(acot_deriv).
                        series_func(acot_series).
			set_name("arccot", "\\arccot"));


} // namespace GiNaC
