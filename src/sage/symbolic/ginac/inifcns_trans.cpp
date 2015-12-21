/** @file inifcns_trans.cpp
 *
 *  Implementation of transcendental (and trigonometric and hyperbolic)
 *  functions. */

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


/* In Sage all the arc trig functions are printed with "arc" instead
   of "a" at the beginning.   This is for consistency with other
   computer algebra systems.   These print methods are registered
   below with each of the corresponding inverse trig function. */


//////////
// exponential function
//////////

static ex exp_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return exp(ex_to<numeric>(x));
	
	return exp(x).hold();
}

static ex exp_eval(const ex & x)
{
	// exp(0) -> 1
	if (x.is_zero()) {
		return _ex1;
	}

	// handle infinity
	// This needs to be before the other tests below, since multiplying
	// infinity with other values throws runtime_errors
	// exp(oo) -> oo
	// exp(-oo) -> 0
	// exp(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity))
			return Infinity;
		if (x.is_equal(NegInfinity))
			return _ex0;
		// x is UnsignedInfinity
		throw (std::runtime_error("exp_eval(): exp^(unsigned_infinity) encountered"));
	}

        bool has_pi = false, has_py = false;
        std::function<bool (const ex&)> has_pi_and_py =
        [&](const ex & the_ex) -> bool
        {
                if (is_exactly_a<constant>(the_ex)
                        and ex_to<constant>(the_ex) == Pi)
                        has_pi = true;
                if (is_a_python_object(the_ex))
                        has_py = true;
                if (has_pi and has_py)
                        return true;
                for (size_t i=0; i<the_ex.nops(); ++i)
                        if (has_pi_and_py(the_ex.op(i)))
                                return true;
                return false;
        };

        if (has_pi_and_py(x)) { // true if x contains Pi and Python objects
                                // like I. The following
                // exp(n*Pi*I/2) -> {+1|+I|-1|-I}
                ex TwoExOverPiI;
                // Arithmetic with I can result in Python exceptions
                // We just ignore the error in this case and skip this step
                try {
                        TwoExOverPiI=(_ex2*x)/(Pi*I);
                } catch (...) {
                }
                if (PyErr_Occurred()) {
                        PyErr_Clear();
                } else if (is_exactly_a<numeric>(TwoExOverPiI)
                        and TwoExOverPiI.info(info_flags::integer)) {
                        const numeric z = mod(ex_to<numeric>(TwoExOverPiI),*_num4_p);
                        if (z.is_equal(*_num0_p))
                                return _ex1;
                        if (z.is_equal(*_num1_p))
                                return ex(I);
                        if (z.is_equal(*_num2_p))
                                return _ex_1;
                        if (z.is_equal(*_num3_p))
                                return ex(-I);
                }
        }

	// exp(log(x)) -> x
	if (is_ex_the_function(x, log))
		return x.op(0);
	
	// exp(float) -> float
	if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
		return exp(ex_to<numeric>(x));
	
	return exp(x).hold();
}

static ex exp_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx exp(x) -> exp(x)
	return exp(x);
}

static ex exp_real_part(const ex & x)
{
	return exp(GiNaC::real_part(x))*cos(GiNaC::imag_part(x));
}

static ex exp_imag_part(const ex & x)
{
	return exp(GiNaC::real_part(x))*sin(GiNaC::imag_part(x));
}

static ex exp_power(const ex & arg, const ex & p)
{
	// If you change this function make sure to adjust
	// the simplification code in mul::eval accordingly
	if (is_exactly_a<numeric>(p) && ex_to<numeric>(p).is_integer())
		return exp(p*arg);
	else
		return power(exp(arg), p).hold();
}

static void exp_print(const ex & arg, const print_context & c,
		bool latex=false)
{
	c.s << "e";
	if (!arg.is_equal(*_num1_p)) {
		c.s<<"^";
		std::stringstream tstream;
		std::unique_ptr<print_context> tcontext_p;
		if (latex) {
			tcontext_p.reset(new print_latex(tstream, c.options));
		} else {
			tcontext_p.reset(new print_dflt(tstream, c.options));
		}
		arg.print(*tcontext_p);
		std::string argstr = tstream.str();
		bool parenthesis = ((argstr.find(' ') != std::string::npos)||
				(argstr.find('+') != std::string::npos) ||
				(argstr.find('-') != std::string::npos) ||
				(argstr.find('/') != std::string::npos) ||
				(argstr.find('*') != std::string::npos) ||
				(argstr.find('^') != std::string::npos));
		if (latex) {
			c.s << '{';
			if (parenthesis)
				c.s << "\\left(";
		} else if (parenthesis)
			c.s << '(';

		c.s << argstr;
		if (latex) {
			if (parenthesis)
				c.s << "\\right)";
			c.s << '}';
		} else if (parenthesis)
			c.s << ')';
	}
}

static void exp_print_dflt(const ex & arg, const print_context & c)
{
	exp_print(arg, c, false);
}

static void exp_print_latex(const ex & arg, const print_context & c)
{
	exp_print(arg, c, true);
}

static ex exp_conjugate(const ex & x)
{
	// conjugate(exp(x))==exp(conjugate(x))
	return exp(x.conjugate());
}

REGISTER_FUNCTION(exp, eval_func(exp_eval).
                       evalf_func(exp_evalf).
                       derivative_func(exp_deriv).
                       real_part_func(exp_real_part).
                       imag_part_func(exp_imag_part).
                       power_func(exp_power).
                       conjugate_func(exp_conjugate).
                       print_func<print_dflt>(exp_print_dflt).
                       print_func<print_latex>(exp_print_latex));

//////////
// natural logarithm
//////////

static ex log_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return log(ex_to<numeric>(x));
	
	return log(x).hold();
}

static ex log_eval(const ex & x)
{
	if (x.info(info_flags::numeric)) {
		// log(float) -> float
		if (!x.info(info_flags::crational))
			return log(ex_to<numeric>(x));

		if (x.is_zero())         // log(0) -> infinity
			//throw(pole_error("log_eval(): log(0)",0));
			return NegInfinity;
		if (not x.info(info_flags::inexact) and x.info(info_flags::negative))
			return (log(-x)+I*Pi);
		if (x.is_equal(_ex1))  // log(1) -> 0
			return _ex0;
		if (x.is_equal(I))       // log(I) -> Pi*I/2
			return (Pi*I*_ex1_2);
		if (x.is_equal(-I))      // log(-I) -> -Pi*I/2
			return (Pi*I*_ex_1_2);

	}

	// log(exp(t)) -> t (if -Pi < t.imag() <= Pi):
	if (is_ex_the_function(x, exp)) {
		const ex &t = x.op(0);
		if (t.info(info_flags::real))
			return t;
	}
	
	// log(oo) -> oo
	if (x.info(info_flags::infinity)) {
		return Infinity;
	}

	return log(x).hold();
}

static ex log_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx log(x) -> 1/x
	return power(x, _ex_1);
}

static ex log_series(const ex &arg,
                     const relational &rel,
                     int order,
                     unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	ex arg_pt;
	bool must_expand_arg = false;
	// maybe substitution of rel into arg fails because of a pole
	try {
		arg_pt = arg.subs(rel, subs_options::no_pattern);
	} catch (pole_error) {
		must_expand_arg = true;
	}
	// or we are at the branch point anyways
	if (arg_pt.is_zero())
		must_expand_arg = true;
	
	if (arg.diff(ex_to<symbol>(rel.lhs())).is_zero()) {
		throw do_taylor();
	}

	if (must_expand_arg) {
		// method:
		// This is the branch point: Series expand the argument first, then
		// trivially factorize it to isolate that part which has constant
		// leading coefficient in this fashion:
		//   x^n + x^(n+1) +...+ Order(x^(n+m))  ->  x^n * (1 + x +...+ Order(x^m)).
		// Return a plain n*log(x) for the x^n part and series expand the
		// other part.  Add them together and reexpand again in order to have
		// one unnested pseries object.  All this also works for negative n.
		pseries argser;          // series expansion of log's argument
		unsigned extra_ord = 0;  // extra expansion order
		do {
			// oops, the argument expanded to a pure Order(x^something)...
			argser = ex_to<pseries>(arg.series(rel, order+extra_ord, options));
			++extra_ord;
		} while (!argser.is_terminating() && argser.nops()==1);

		const symbol &s = ex_to<symbol>(rel.lhs());
		const ex &point = rel.rhs();
		const int n = argser.ldegree(s);
		epvector seq;
		// construct what we carelessly called the n*log(x) term above
		const ex coeff = argser.coeff(s, n);
		// expand the log, but only if coeff is real and > 0, since otherwise
		// it would make the branch cut run into the wrong direction
		if (coeff.info(info_flags::positive))
			seq.push_back(expair(n*log(s-point)+log(coeff), _ex0));
		else
			seq.push_back(expair(log(coeff*pow(s-point, n)), _ex0));

		if (!argser.is_terminating() || argser.nops()!=1) {
			// in this case n more (or less) terms are needed
			// (sadly, to generate them, we have to start from the beginning)
			if (n == 0 && coeff == 1) {
				epvector epv;
				ex acc = (new pseries(rel, epv))->setflag(status_flags::dynallocated);
				epv.reserve(2);
				epv.push_back(expair(-1, _ex0));
				epv.push_back(expair(Order(_ex1), order));
				ex rest = pseries(rel, epv).add_series(argser);
				for (int i = order-1; i>0; --i) {
					epvector cterm;
					cterm.reserve(1);
					cterm.push_back(expair(i%2 ? _ex1/i : _ex_1/i, _ex0));
					acc = pseries(rel, cterm).add_series(ex_to<pseries>(acc));
					acc = (ex_to<pseries>(rest)).mul_series(ex_to<pseries>(acc));
				}
				return acc;
			}
			const ex newarg = ex_to<pseries>((arg/coeff).series(rel, order+n, options)).shift_exponents(-n).convert_to_poly(true);
			return pseries(rel, seq).add_series(ex_to<pseries>(log(newarg).series(rel, order, options)));
		} else  // it was a monomial
			return pseries(rel, seq);
	}
	if (!(options & series_options::suppress_branchcut) &&
	     arg_pt.info(info_flags::negative)) {
		// method:
		// This is the branch cut: assemble the primitive series manually and
		// then add the corresponding complex step function.
		const symbol &s = ex_to<symbol>(rel.lhs());
		const ex &point = rel.rhs();
		const symbol foo;
		const ex replarg = series(log(arg), s==foo, order).subs(foo==point, subs_options::no_pattern);
		epvector seq;
		seq.push_back(expair(-I*csgn(arg*I)*Pi, _ex0));
		seq.push_back(expair(Order(_ex1), order));
		return series(replarg - I*Pi + pseries(rel, seq), rel, order);
	}
	throw do_taylor();  // caught by function::series()
}

static ex log_real_part(const ex & x)
{
	if (x.info(info_flags::positive))
		return log(x).hold();
	return log(abs(x));
}

static ex log_imag_part(const ex & x)
{
	if (x.info(info_flags::positive))
		return 0;
	return atan2(GiNaC::imag_part(x), GiNaC::real_part(x));
}

static ex log_conjugate(const ex & x)
{
	// conjugate(log(x))==log(conjugate(x)) unless on the branch cut which
	// runs along the negative real axis.
	if (x.info(info_flags::positive)) {
		return log(x);
	}
	if (is_exactly_a<numeric>(x) &&
	    !x.imag_part().is_zero()) {
		return log(x.conjugate());
	}
	return conjugate_function(log(x)).hold();
}

REGISTER_FUNCTION(log, eval_func(log_eval).
                       evalf_func(log_evalf).
                       derivative_func(log_deriv).
                       series_func(log_series).
                       real_part_func(log_real_part).
                       imag_part_func(log_imag_part).
                       conjugate_func(log_conjugate).
                       latex_name("\\log"));

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
	if (TwentyforExOverPi.info(info_flags::integer)) {
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
	if (FortyeightExOverPi.info(info_flags::integer)) {
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
	ex a = GiNaC::real_part(x);
	ex b = GiNaC::imag_part(x);
	return tan(a)/(1+power(tan(a),2)*power(tan(b),2));
}

static ex tan_imag_part(const ex & x)
{
	ex a = GiNaC::real_part(x);
	ex b = GiNaC::imag_part(x);
	return tanh(b)/(1+power(tan(a),2)*power(tan(b),2));
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
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity))
			return _ex1_2*Pi;
		if (x.is_equal(NegInfinity))
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
	if (x.info(info_flags::numeric)) {

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
	
	if ((x/Pi).info(info_flags::numeric) &&
		ex_to<numeric>(x/Pi).real().is_zero())  // sinh(I*x) -> I*sin(x)
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
	if (x.info(info_flags::numeric)) {

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

	if ((x/Pi).info(info_flags::numeric) &&
		ex_to<numeric>(x/Pi).real().is_zero())  // cosh(I*x) -> cos(x)
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
	if (x.info(info_flags::numeric)) {

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
		
	if ((x/Pi).info(info_flags::numeric) &&
		ex_to<numeric>(x/Pi).real().is_zero())  // tanh(I*x) -> I*tan(x);
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

static ex tanh_real_part(const ex & x)
{
	ex a = GiNaC::real_part(x);
	ex b = GiNaC::imag_part(x);
	return tanh(a)/(1+power(tanh(a),2)*power(tan(b),2));
}

static ex tanh_imag_part(const ex & x)
{
	ex a = GiNaC::real_part(x);
	ex b = GiNaC::imag_part(x);
	return tan(b)/(1+power(tanh(a),2)*power(tan(b),2));
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
	if (x.info(info_flags::numeric)) {

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
	if (x.info(info_flags::numeric)) {

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
	if (x.info(info_flags::numeric)) {

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


} // namespace GiNaC
