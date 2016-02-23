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

	ex x_red;
        if (has_pi_and_py(x)) { // true if x contains Pi and Python objects
                                // like I. To make this check should be faster
                                // than the following

		ex coef_pi = x.coeff(Pi).expand();
		ex rem = _ex0;
		if (is_exactly_a<add>(coef_pi)) {
			for (size_t i=0; i < coef_pi.nops(); i++) {
				if ((coef_pi.op(i) / (_ex2 * I)).info(info_flags::integer))
					rem += Pi * coef_pi.op(i);
			}
		}
		else if ((coef_pi / (_ex2 * I)).info(info_flags::integer))
			rem = Pi * coef_pi;
		x_red = x - rem;


		ex res1 = sin(x_red/I);
		ex res2 = cos(x_red/I);
		if (((not is_exactly_a<function>(res1))
		or (not is_ex_the_function(res1, sin)))
		and ((not is_exactly_a<function>(res2))
		or (not is_ex_the_function(res2, cos))))
			return I*res1 + res2;
        }
	else
		x_red = x;

	// exp(log(x)) -> x
	if (is_ex_the_function(x_red, log))
		return x_red.op(0);
	
	// exp(float) -> float
	if (is_exactly_a<numeric>(x_red) && !x_red.info(info_flags::crational))
		return exp(ex_to<numeric>(x_red));
	
	return exp(x_red).hold();
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
	if (is_exactly_a<numeric>(x)) {
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

} // namespace GiNaC
