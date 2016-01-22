/** @file inifcns_trans.cpp
 *
 *  Implementation of transcendental (and hyperbolic) functions. */

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
	if (is_exactly_a<infinity>(x)) {
	        infinity xinf = ex_to<infinity>(x);
		if (xinf.is_plus_infinity())
			return Infinity;
		if (xinf.is_minus_infinity())
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
                                // like I. To make this check should be faster
                                // than the following
                ex res1 = sin(x/I);
                ex res2 = cos(x/I);
                if (((not is_exactly_a<function>(res1))
                or (not is_ex_the_function(res1, sin)))
                and ((not is_exactly_a<function>(res2))
                or (not is_ex_the_function(res2, cos))))
                        return I*res1 + res2;
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
