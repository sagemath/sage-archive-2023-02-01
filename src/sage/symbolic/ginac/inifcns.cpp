/** @file inifcns.cpp
 *
 *  Implementation of GiNaC's initially known functions. */

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

#include <vector>
#include <stdexcept>

#include "inifcns.h"
#include "ex.h"
#include "constant.h"
#include "lst.h"
#include "matrix.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "pseries.h"
#include "symbol.h"
#include "symmetry.h"
#include "utils.h"

namespace GiNaC {

//////////
// complex conjugate
//////////

static ex conjugate_evalf(const ex & arg, int prec)
{
	if (is_exactly_a<numeric>(arg)) {
		return ex_to<numeric>(arg).conjugate();
	}
	return conjugate_function(arg).hold();
}

static ex conjugate_eval(const ex & arg)
{
	return arg.conjugate();
}

static void conjugate_print_latex(const ex & arg, const print_context & c)
{
	c.s << "\\bar{"; arg.print(c); c.s << "}";
}

static ex conjugate_conjugate(const ex & arg)
{
	return arg;
}

static ex conjugate_real_part(const ex & arg)
{
	return arg.real_part();
}

static ex conjugate_imag_part(const ex & arg)
{
	return -arg.imag_part();
}

REGISTER_FUNCTION(conjugate_function, eval_func(conjugate_eval).
                                      evalf_func(conjugate_evalf).
                                      print_func<print_latex>(conjugate_print_latex).
                                      conjugate_func(conjugate_conjugate).
                                      real_part_func(conjugate_real_part).
                                      imag_part_func(conjugate_imag_part).
                                      set_name("conjugate","conjugate"));

//////////
// real part
//////////

static ex real_part_evalf(const ex & arg, int prec)
{
	if (is_exactly_a<numeric>(arg)) {
		return ex_to<numeric>(arg).real();
	}
	return real_part_function(arg).hold();
}

static ex real_part_eval(const ex & arg)
{
	return arg.real_part();
}

static void real_part_print_latex(const ex & arg, const print_context & c)
{
	c.s << "\\Re"; arg.print(c); c.s << "";
}

static ex real_part_conjugate(const ex & arg)
{
	return real_part_function(arg).hold();
}

static ex real_part_real_part(const ex & arg)
{
	return real_part_function(arg).hold();
}

static ex real_part_imag_part(const ex & arg)
{
	return 0;
}

REGISTER_FUNCTION(real_part_function, eval_func(real_part_eval).
                                      evalf_func(real_part_evalf).
                                      print_func<print_latex>(real_part_print_latex).
                                      conjugate_func(real_part_conjugate).
                                      real_part_func(real_part_real_part).
                                      imag_part_func(real_part_imag_part).
                                      set_name("real_part","real_part"));

//////////
// imag part
//////////

static ex imag_part_evalf(const ex & arg, int prec)
{
	if (is_exactly_a<numeric>(arg)) {
		return ex_to<numeric>(arg).imag();
	}
	return imag_part_function(arg).hold();
}

static ex imag_part_eval(const ex & arg)
{
	return arg.imag_part();
}

static void imag_part_print_latex(const ex & arg, const print_context & c)
{
	c.s << "\\Im"; arg.print(c); c.s << "";
}

static ex imag_part_conjugate(const ex & arg)
{
	return imag_part_function(arg).hold();
}

static ex imag_part_real_part(const ex & arg)
{
	return imag_part_function(arg).hold();
}

static ex imag_part_imag_part(const ex & arg)
{
	return 0;
}

REGISTER_FUNCTION(imag_part_function, eval_func(imag_part_eval).
                                      evalf_func(imag_part_evalf).
                                      print_func<print_latex>(imag_part_print_latex).
                                      conjugate_func(imag_part_conjugate).
                                      real_part_func(imag_part_real_part).
                                      imag_part_func(imag_part_imag_part).
                                      set_name("imag_part","imag_part"));

//////////
// absolute value
//////////

static ex abs_evalf(const ex & arg, int prec)
{
	if (is_exactly_a<numeric>(arg))
		return abs(ex_to<numeric>(arg));
	
	return abs(arg).hold();
}

static ex abs_eval(const ex & arg)
{
	if (is_exactly_a<numeric>(arg))
		return abs(ex_to<numeric>(arg));

	if (arg.info(info_flags::nonnegative))
		return arg;

	if (is_ex_the_function(arg, abs))
		return arg;

	return abs(arg).hold();
}

static void abs_print_latex(const ex & arg, const print_context & c)
{
	c.s << "{\\left|"; arg.print(c); c.s << "\\right|}";
}

static void abs_print_csrc_float(const ex & arg, const print_context & c)
{
	c.s << "fabs("; arg.print(c); c.s << ")";
}

static ex abs_conjugate(const ex & arg)
{
	return abs(arg);
}

static ex abs_real_part(const ex & arg)
{
	return abs(arg).hold();
}

static ex abs_imag_part(const ex& arg)
{
	return 0;
}

static ex abs_power(const ex & arg, const ex & exp)
{
	if (arg.is_equal(arg.conjugate()) && is_a<numeric>(exp) && ex_to<numeric>(exp).is_even())
		return power(arg, exp);
	else
		return power(abs(arg), exp).hold();
}

REGISTER_FUNCTION(abs, eval_func(abs_eval).
                       evalf_func(abs_evalf).
                       print_func<print_latex>(abs_print_latex).
                       print_func<print_csrc_float>(abs_print_csrc_float).
                       print_func<print_csrc_double>(abs_print_csrc_float).
                       conjugate_func(abs_conjugate).
                       real_part_func(abs_real_part).
                       imag_part_func(abs_imag_part).
                       power_func(abs_power));

//////////
// Step function
//////////

static ex step_evalf(const ex & arg, int prec)
{
	if (is_exactly_a<numeric>(arg))
		return step(ex_to<numeric>(arg));
	
	return step(arg).hold();
}

static ex step_eval(const ex & arg)
{
	if (is_exactly_a<numeric>(arg))
		return step(ex_to<numeric>(arg));
	
	else if (is_exactly_a<mul>(arg) &&
	         is_exactly_a<numeric>(arg.op(arg.nops()-1))) {
		numeric oc = ex_to<numeric>(arg.op(arg.nops()-1));
		if (oc.is_real()) {
			if (oc > 0)
				// step(42*x) -> step(x)
				return step(arg/oc).hold();
			else
				// step(-42*x) -> step(-x)
				return step(-arg/oc).hold();
		}
		if (oc.real().is_zero()) {
			if (oc.imag() > 0)
				// step(42*I*x) -> step(I*x)
				return step(I*arg/oc).hold();
			else
				// step(-42*I*x) -> step(-I*x)
				return step(-I*arg/oc).hold();
		}
	}
	
	return step(arg).hold();
}

static ex step_series(const ex & arg,
                      const relational & rel,
                      int order,
                      unsigned options)
{
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (arg_pt.info(info_flags::numeric)
	    && ex_to<numeric>(arg_pt).real().is_zero()
	    && !(options & series_options::suppress_branchcut))
		throw (std::domain_error("step_series(): on imaginary axis"));
	
	epvector seq;
	seq.push_back(expair(step(arg_pt), _ex0));
	return pseries(rel,seq);
}

static ex step_conjugate(const ex& arg)
{
	return step(arg).hold();
}

static ex step_real_part(const ex& arg)
{
	return step(arg).hold();
}

static ex step_imag_part(const ex& arg)
{
	return 0;
}

REGISTER_FUNCTION(step, eval_func(step_eval).
                        evalf_func(step_evalf).
                        series_func(step_series).
                        conjugate_func(step_conjugate).
								real_part_func(step_real_part).
								imag_part_func(step_imag_part));

//////////
// Complex sign
//////////

static ex csgn_evalf(const ex & arg, int prec)
{
	if (is_exactly_a<numeric>(arg))
		return csgn(ex_to<numeric>(arg));
	
	return csgn(arg).hold();
}

static ex csgn_eval(const ex & arg)
{
	if (is_exactly_a<numeric>(arg))
		return csgn(ex_to<numeric>(arg));
	
	else if (is_exactly_a<mul>(arg) &&
	         is_exactly_a<numeric>(arg.op(arg.nops()-1))) {
		numeric oc = ex_to<numeric>(arg.op(arg.nops()-1));
		if (oc.is_real()) {
			if (oc > 0)
				// csgn(42*x) -> csgn(x)
				return csgn(arg/oc).hold();
			else
				// csgn(-42*x) -> -csgn(x)
				return -csgn(arg/oc).hold();
		}
		if (oc.real().is_zero()) {
			if (oc.imag() > 0)
				// csgn(42*I*x) -> csgn(I*x)
				return csgn(I*arg/oc).hold();
			else
				// csgn(-42*I*x) -> -csgn(I*x)
				return -csgn(I*arg/oc).hold();
		}
	}
	
	return csgn(arg).hold();
}

static ex csgn_series(const ex & arg,
                      const relational & rel,
                      int order,
                      unsigned options)
{
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (arg_pt.info(info_flags::numeric)
	    && ex_to<numeric>(arg_pt).real().is_zero()
	    && !(options & series_options::suppress_branchcut))
		throw (std::domain_error("csgn_series(): on imaginary axis"));
	
	epvector seq;
	seq.push_back(expair(csgn(arg_pt), _ex0));
	return pseries(rel,seq);
}

static ex csgn_conjugate(const ex& arg)
{
	return csgn(arg).hold();
}

static ex csgn_real_part(const ex& arg)
{
	return csgn(arg).hold();
}

static ex csgn_imag_part(const ex& arg)
{
	return 0;
}

static ex csgn_power(const ex & arg, const ex & exp)
{
	if (is_a<numeric>(exp) && exp.info(info_flags::positive) && ex_to<numeric>(exp).is_integer()) {
		if (ex_to<numeric>(exp).is_odd())
			return csgn(arg);
		else
			return power(csgn(arg), _ex2).hold();
	} else
		return power(csgn(arg), exp).hold();
}


REGISTER_FUNCTION(csgn, eval_func(csgn_eval).
                        evalf_func(csgn_evalf).
                        series_func(csgn_series).
                        conjugate_func(csgn_conjugate).
                        real_part_func(csgn_real_part).
                        imag_part_func(csgn_imag_part).
                        power_func(csgn_power));


//////////
// Eta function: eta(x,y) == log(x*y) - log(x) - log(y).
// This function is closely related to the unwinding number K, sometimes found
// in modern literature: K(z) == (z-log(exp(z)))/(2*Pi*I).
//////////

static ex eta_evalf(const ex &x, const ex &y, int prec)
{
	// It seems like we basically have to replicate the eval function here,
	// since the expression might not be fully evaluated yet.
	if (x.info(info_flags::positive) || y.info(info_flags::positive))
		return _ex0;

	if (x.info(info_flags::numeric) &&	y.info(info_flags::numeric)) {
		const numeric nx = ex_to<numeric>(x);
		const numeric ny = ex_to<numeric>(y);
		const numeric nxy = ex_to<numeric>(x*y);
		int cut = 0;
		if (nx.is_real() && nx.is_negative())
			cut -= 4;
		if (ny.is_real() && ny.is_negative())
			cut -= 4;
		if (nxy.is_real() && nxy.is_negative())
			cut += 4;
		return evalf(I/4*Pi, 0, prec)*((csgn(-imag(nx))+1)*(csgn(-imag(ny))+1)*(csgn(imag(nxy))+1)-
		                      (csgn(imag(nx))+1)*(csgn(imag(ny))+1)*(csgn(-imag(nxy))+1)+cut);
	}

	return eta(x,y).hold();
}

static ex eta_eval(const ex &x, const ex &y)
{
	// trivial:  eta(x,c) -> 0  if c is real and positive
	if (x.info(info_flags::positive) || y.info(info_flags::positive))
		return _ex0;

	if (x.info(info_flags::numeric) &&	y.info(info_flags::numeric)) {
		// don't call eta_evalf here because it would call Pi.evalf()!
		const numeric nx = ex_to<numeric>(x);
		const numeric ny = ex_to<numeric>(y);
		const numeric nxy = ex_to<numeric>(x*y);
		int cut = 0;
		if (nx.is_real() && nx.is_negative())
			cut -= 4;
		if (ny.is_real() && ny.is_negative())
			cut -= 4;
		if (nxy.is_real() && nxy.is_negative())
			cut += 4;
		return (I/4)*Pi*((csgn(-imag(nx))+1)*(csgn(-imag(ny))+1)*(csgn(imag(nxy))+1)-
		                 (csgn(imag(nx))+1)*(csgn(imag(ny))+1)*(csgn(-imag(nxy))+1)+cut);
	}
	
	return eta(x,y).hold();
}

static ex eta_series(const ex & x, const ex & y,
                     const relational & rel,
                     int order,
                     unsigned options)
{
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	const ex y_pt = y.subs(rel, subs_options::no_pattern);
	if ((x_pt.info(info_flags::numeric) && x_pt.info(info_flags::negative)) ||
	    (y_pt.info(info_flags::numeric) && y_pt.info(info_flags::negative)) ||
	    ((x_pt*y_pt).info(info_flags::numeric) && (x_pt*y_pt).info(info_flags::negative)))
			throw (std::domain_error("eta_series(): on discontinuity"));
	epvector seq;
	seq.push_back(expair(eta(x_pt,y_pt), _ex0));
	return pseries(rel,seq);
}

static ex eta_conjugate(const ex & x, const ex & y)
{
	return -eta(x, y);
}

static ex eta_real_part(const ex & x, const ex & y)
{
	return 0;
}

static ex eta_imag_part(const ex & x, const ex & y)
{
	return -I*eta(x, y).hold();
}

REGISTER_FUNCTION(eta, eval_func(eta_eval).
                       evalf_func(eta_evalf).
                       series_func(eta_series).
                       latex_name("\\eta").
                       set_symmetry(sy_symm(0, 1)).
                       conjugate_func(eta_conjugate).
                       real_part_func(eta_real_part).
                       imag_part_func(eta_imag_part));


//////////
// dilogarithm
//////////

static ex Li2_evalf(const ex & x, int prec)
{
	if (is_exactly_a<numeric>(x))
		return Li2(ex_to<numeric>(x));
	
	return Li2(x).hold();
}

static ex Li2_eval(const ex & x)
{
	if (x.info(info_flags::numeric)) {
		// Li2(0) -> 0
		if (x.is_zero())
			return _ex0;
		// Li2(1) -> Pi^2/6
		if (x.is_equal(_ex1))
			return power(Pi,_ex2)/_ex6;
		// Li2(1/2) -> Pi^2/12 - log(2)^2/2
		if (x.is_equal(_ex1_2))
			return power(Pi,_ex2)/_ex12 + power(log(_ex2),_ex2)*_ex_1_2;
		// Li2(-1) -> -Pi^2/12
		if (x.is_equal(_ex_1))
			return -power(Pi,_ex2)/_ex12;
		// Li2(I) -> -Pi^2/48+Catalan*I
		if (x.is_equal(I))
			return power(Pi,_ex2)/_ex_48 + Catalan*I;
		// Li2(-I) -> -Pi^2/48-Catalan*I
		if (x.is_equal(-I))
			return power(Pi,_ex2)/_ex_48 - Catalan*I;
		// Li2(float)
		/* delay numeric evaluation to evalf,
		 * where precision can be specified
		if (!x.info(info_flags::crational))
			return Li2(ex_to<numeric>(x));
		*/
	}
	
	return Li2(x).hold();
}

static ex Li2_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx Li2(x) -> -log(1-x)/x
	return -log(_ex1-x)/x;
}

static ex Li2_series(const ex &x, const relational &rel, int order, unsigned options)
{
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (x_pt.info(info_flags::numeric)) {
		// First special case: x==0 (derivatives have poles)
		if (x_pt.is_zero()) {
			// method:
			// The problem is that in d/dx Li2(x==0) == -log(1-x)/x we cannot 
			// simply substitute x==0.  The limit, however, exists: it is 1.
			// We also know all higher derivatives' limits:
			// (d/dx)^n Li2(x) == n!/n^2.
			// So the primitive series expansion is
			// Li2(x==0) == x + x^2/4 + x^3/9 + ...
			// and so on.
			// We first construct such a primitive series expansion manually in
			// a dummy symbol s and then insert the argument's series expansion
			// for s.  Reexpanding the resulting series returns the desired
			// result.
			const symbol s;
			ex ser;
			// manually construct the primitive expansion
			for (int i=1; i<order; ++i)
				ser += pow(s,i) / pow(numeric(i), *_num2_p);
			// substitute the argument's series expansion
			ser = ser.subs(s==x.series(rel, order), subs_options::no_pattern);
			// maybe that was terminating, so add a proper order term
			epvector nseq;
			nseq.push_back(expair(Order(_ex1), order));
			ser += pseries(rel, nseq);
			// reexpanding it will collapse the series again
			return ser.series(rel, order);
			// NB: Of course, this still does not allow us to compute anything
			// like sin(Li2(x)).series(x==0,2), since then this code here is
			// not reached and the derivative of sin(Li2(x)) doesn't allow the
			// substitution x==0.  Probably limits *are* needed for the general
			// cases.  In case L'Hospital's rule is implemented for limits and
			// basic::series() takes care of this, this whole block is probably
			// obsolete!
		}
		// second special case: x==1 (branch point)
		if (x_pt.is_equal(_ex1)) {
			// method:
			// construct series manually in a dummy symbol s
			const symbol s;
			ex ser = zeta(_ex2);
			// manually construct the primitive expansion
			for (int i=1; i<order; ++i)
				ser += pow(1-s,i) * (numeric(1,i)*(I*Pi+log(s-1)) - numeric(1,i*i));
			// substitute the argument's series expansion
			ser = ser.subs(s==x.series(rel, order), subs_options::no_pattern);
			// maybe that was terminating, so add a proper order term
			epvector nseq;
			nseq.push_back(expair(Order(_ex1), order));
			ser += pseries(rel, nseq);
			// reexpanding it will collapse the series again
			return ser.series(rel, order);
		}
		// third special case: x real, >=1 (branch cut)
		if (!(options & series_options::suppress_branchcut) &&
			ex_to<numeric>(x_pt).is_real() && ex_to<numeric>(x_pt)>1) {
			// method:
			// This is the branch cut: assemble the primitive series manually
			// and then add the corresponding complex step function.
			const symbol &s = ex_to<symbol>(rel.lhs());
			const ex point = rel.rhs();
			const symbol foo;
			epvector seq;
			// zeroth order term:
			seq.push_back(expair(Li2(x_pt), _ex0));
			// compute the intermediate terms:
			ex replarg = series(Li2(x), s==foo, order);
			for (size_t i=1; i<replarg.nops()-1; ++i)
				seq.push_back(expair((replarg.op(i)/power(s-foo,i)).series(foo==point,1,options).op(0).subs(foo==s, subs_options::no_pattern),i));
			// append an order term:
			seq.push_back(expair(Order(_ex1), replarg.nops()-1));
			return pseries(rel, seq);
		}
	}
	// all other cases should be safe, by now:
	throw do_taylor();  // caught by function::series()
}

REGISTER_FUNCTION(Li2, eval_func(Li2_eval).
                       evalf_func(Li2_evalf).
                       derivative_func(Li2_deriv).
                       series_func(Li2_series).
                       latex_name("\\mbox{Li}_2"));

//////////
// trilogarithm
//////////

static ex Li3_eval(const ex & x)
{
	if (x.is_zero())
		return x;
	return Li3(x).hold();
}

REGISTER_FUNCTION(Li3, eval_func(Li3_eval).
                       latex_name("\\mbox{Li}_3"));

//////////
// Derivatives of Riemann's Zeta-function  zetaderiv(0,x)==zeta(x)
//////////

static ex zetaderiv_eval(const ex & n, const ex & x)
{
	if (n.info(info_flags::numeric)) {
		// zetaderiv(0,x) -> zeta(x)
		if (n.is_zero())
			return zeta(x);
	}
	
	return zetaderiv(n, x).hold();
}

static ex zetaderiv_deriv(const ex & n, const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param<2);
	
	if (deriv_param==0) {
		// d/dn zeta(n,x)
		throw(std::logic_error("cannot diff zetaderiv(n,x) with respect to n"));
	}
	// d/dx psi(n,x)
	return zetaderiv(n+1,x);
}

REGISTER_FUNCTION(zetaderiv, eval_func(zetaderiv_eval).
	                       	 derivative_func(zetaderiv_deriv).
  	                         latex_name("\\zeta^\\prime"));

//////////
// factorial
//////////

static ex factorial_evalf(const ex & x, int prec)
{
	return factorial(x).hold();
}

static ex factorial_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x))
		return factorial(ex_to<numeric>(x));
	else
		return factorial(x).hold();
}

static void factorial_print_dflt_latex(const ex & x, const print_context & c)
{
	if (is_exactly_a<symbol>(x) ||
	    is_exactly_a<constant>(x) ||
		is_exactly_a<function>(x)) {
		x.print(c); c.s << "!";
	} else {
		c.s << "\\left("; x.print(c); c.s << "\\right)!";
	}
}

static ex factorial_conjugate(const ex & x)
{
	return factorial(x).hold();
}

static ex factorial_real_part(const ex & x)
{
	return factorial(x).hold();
}

static ex factorial_imag_part(const ex & x)
{
	return 0;
}

REGISTER_FUNCTION(factorial, eval_func(factorial_eval).
                             evalf_func(factorial_evalf).
                           //print_func<print_dflt>(factorial_print_dflt_latex).
                             print_func<print_latex>(factorial_print_dflt_latex).
                             conjugate_func(factorial_conjugate).
                             real_part_func(factorial_real_part).
                             imag_part_func(factorial_imag_part));

//////////
// binomial
//////////

static ex binomial_evalf(const ex & x, const ex & y, int prec)
{
	return binomial(x, y).hold();
}

static ex binomial_sym(const ex & x, const numeric & y)
{
	if (y.is_integer()) {
		if (y.is_nonneg_integer()) {
			const unsigned N = y.to_int();
			if (N == 0) return _ex0;
			if (N == 1) return x;
			ex t = x.expand();
			for (unsigned i = 2; i <= N; ++i)
				t = (t * (x + i - y - 1)).expand() / i;
			return t;
		} else
			return _ex0;
	}

	return binomial(x, y).hold();
}

static ex binomial_eval(const ex & x, const ex &y)
{
	if (is_exactly_a<numeric>(y)) {
		if (is_exactly_a<numeric>(x) && ex_to<numeric>(x).is_integer())
			return binomial(ex_to<numeric>(x), ex_to<numeric>(y));
		else
			return binomial_sym(x, ex_to<numeric>(y));
	} else
		return binomial(x, y).hold();
}

// At the moment the numeric evaluation of a binomail function always
// gives a real number, but if this would be implemented using the gamma
// function, also complex conjugation should be changed (or rather, deleted).
static ex binomial_conjugate(const ex & x, const ex & y)
{
	return binomial(x,y).hold();
}

static ex binomial_real_part(const ex & x, const ex & y)
{
	return binomial(x,y).hold();
}

static ex binomial_imag_part(const ex & x, const ex & y)
{
	return 0;
}

static void binomial_print_latex(const ex & x, const ex & y,
		const print_context & c)
{
	c.s<<"{";
	x.print(c);
	c.s<<" \\choose ";
	y.print(c);
	c.s<<"}";
}


REGISTER_FUNCTION(binomial, eval_func(binomial_eval).
                            evalf_func(binomial_evalf).
                            conjugate_func(binomial_conjugate).
                            real_part_func(binomial_real_part).
                            print_func<print_latex>(binomial_print_latex).
                            imag_part_func(binomial_imag_part));

//////////
// Order term function (for truncated power series)
//////////

static ex Order_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {
		// O(c) -> O(1) or 0
		if (!x.is_zero())
			return Order(_ex1).hold();
		else
			return _ex0;
	} else if (is_exactly_a<mul>(x)) {
		const mul &m = ex_to<mul>(x);
		// O(c*expr) -> O(expr)
		if (is_exactly_a<numeric>(m.op(m.nops() - 1)))
			return Order(x / m.op(m.nops() - 1)).hold();
	}
	return Order(x).hold();
}

static ex Order_series(const ex & x, const relational & r, int order, unsigned options)
{
	// Just wrap the function into a pseries object
	epvector new_seq;
	GINAC_ASSERT(is_a<symbol>(r.lhs()));
	const symbol &s = ex_to<symbol>(r.lhs());
	new_seq.push_back(expair(Order(_ex1), numeric(std::min(x.ldegree(s), order))));
	return pseries(r, new_seq);
}

static ex Order_conjugate(const ex & x)
{
	return Order(x).hold();
}

static ex Order_real_part(const ex & x)
{
	return Order(x).hold();
}

static ex Order_imag_part(const ex & x)
{
	if(x.info(info_flags::real))
		return 0;
	return Order(x).hold();
}

// Differentiation is handled in function::derivative because of its special requirements

REGISTER_FUNCTION(Order, eval_func(Order_eval).
                         series_func(Order_series).
                         latex_name("\\mathcal{O}").
                         conjugate_func(Order_conjugate).
                         real_part_func(Order_real_part).
                         imag_part_func(Order_imag_part));

//////////
// Solve linear system
//////////

ex lsolve(const ex &eqns, const ex &symbols, unsigned options)
{
	// solve a system of linear equations
	if (eqns.info(info_flags::relation_equal)) {
		if (!symbols.info(info_flags::symbol))
			throw(std::invalid_argument("lsolve(): 2nd argument must be a symbol"));
		const ex sol = lsolve(lst(eqns),lst(symbols));
		
		GINAC_ASSERT(sol.nops()==1);
		GINAC_ASSERT(is_exactly_a<relational>(sol.op(0)));
		
		return sol.op(0).op(1); // return rhs of first solution
	}
	
	// syntax checks
	if (!eqns.info(info_flags::list)) {
		throw(std::invalid_argument("lsolve(): 1st argument must be a list or an equation"));
	}
	for (size_t i=0; i<eqns.nops(); i++) {
		if (!eqns.op(i).info(info_flags::relation_equal)) {
			throw(std::invalid_argument("lsolve(): 1st argument must be a list of equations"));
		}
	}
	if (!symbols.info(info_flags::list)) {
		throw(std::invalid_argument("lsolve(): 2nd argument must be a list or a symbol"));
	}
	for (size_t i=0; i<symbols.nops(); i++) {
		if (!symbols.op(i).info(info_flags::symbol)) {
			throw(std::invalid_argument("lsolve(): 2nd argument must be a list of symbols"));
		}
	}
	
	// build matrix from equation system
	matrix sys(eqns.nops(),symbols.nops());
	matrix rhs(eqns.nops(),1);
	matrix vars(symbols.nops(),1);
	
	for (size_t r=0; r<eqns.nops(); r++) {
		const ex eq = eqns.op(r).op(0)-eqns.op(r).op(1); // lhs-rhs==0
		ex linpart = eq;
		for (size_t c=0; c<symbols.nops(); c++) {
			const ex co = eq.coeff(ex_to<symbol>(symbols.op(c)),1);
			linpart -= co*symbols.op(c);
			sys(r,c) = co;
		}
		linpart = linpart.expand();
		rhs(r,0) = -linpart;
	}
	
	// test if system is linear and fill vars matrix
	for (size_t i=0; i<symbols.nops(); i++) {
		vars(i,0) = symbols.op(i);
		if (sys.has(symbols.op(i)))
			throw(std::logic_error("lsolve: system is not linear"));
		if (rhs.has(symbols.op(i)))
			throw(std::logic_error("lsolve: system is not linear"));
	}
	
	matrix solution;
	try {
		solution = sys.solve(vars,rhs,options);
	} catch (const std::runtime_error & e) {
		// Probably singular matrix or otherwise overdetermined system:
		// It is consistent to return an empty list
		return lst();
	}
	GINAC_ASSERT(solution.cols()==1);
	GINAC_ASSERT(solution.rows()==symbols.nops());
	
	// return list of equations of the form lst(var1==sol1,var2==sol2,...)
	lst sollist;
	for (size_t i=0; i<symbols.nops(); i++)
		sollist.append(symbols.op(i)==solution(i,0));
	
	return sollist;
}

//////////
// Find real root of f(x) numerically
//////////

const numeric
fsolve(const ex& f_in, const symbol& x, const numeric& x1, const numeric& x2,
		int prec)
{
	if (!x1.is_real() || !x2.is_real()) {
		throw std::runtime_error("fsolve(): interval not bounded by real numbers");
	}
	if (x1==x2) {
		throw std::runtime_error("fsolve(): vanishing interval");
	}
	// xx[0] == left interval limit, xx[1] == right interval limit.
	// fx[0] == f(xx[0]), fx[1] == f(xx[1]).
	// We keep the root bracketed: xx[0]<xx[1] and fx[0]*fx[1]<0.
	numeric xx[2] = { x1<x2 ? x1 : x2,
	                  x1<x2 ? x2 : x1 };
	ex f;
	if (is_a<relational>(f_in)) {
		f = f_in.lhs()-f_in.rhs();
	} else {
		f = f_in;
	}
	const ex fx_[2] = { f.subs(x==xx[0]).evalf(0, prec),
	                    f.subs(x==xx[1]).evalf(0, prec) };
	if (!is_a<numeric>(fx_[0]) || !is_a<numeric>(fx_[1])) {
		throw std::runtime_error("fsolve(): function does not evaluate numerically");
	}
	numeric fx[2] = { ex_to<numeric>(fx_[0]),
	                  ex_to<numeric>(fx_[1]) };
	if (!fx[0].is_real() || !fx[1].is_real()) {
		throw std::runtime_error("fsolve(): function evaluates to complex values at interval boundaries");
	}
	if (fx[0]*fx[1]>=0) {
		throw std::runtime_error("fsolve(): function does not change sign at interval boundaries");
	}

	// The Newton-Raphson method has quadratic convergence!  Simply put, it
	// replaces x with x-f(x)/f'(x) at each step.  -f/f' is the delta:
	const ex ff = normal(-f/f.diff(x));
	int side = 0;  // Start at left interval limit.
	numeric xxprev;
	numeric fxprev;
	do {
		xxprev = xx[side];
		fxprev = fx[side];
		xx[side] += ex_to<numeric>(ff.subs(x==xx[side]).evalf(0, prec));
		fx[side] = ex_to<numeric>(f.subs(x==xx[side]).evalf(0, prec));
		if ((side==0 && xx[0]<xxprev) || (side==1 && xx[1]>xxprev) || xx[0]>xx[1]) {
			// Oops, Newton-Raphson method shot out of the interval.
			// Restore, and try again with the other side instead!
			xx[side] = xxprev;
			fx[side] = fxprev;
			side = !side;
			xxprev = xx[side];
			fxprev = fx[side];
			xx[side] += ex_to<numeric>(ff.subs(x==xx[side]).evalf(0, prec));
			fx[side] = ex_to<numeric>(f.subs(x==xx[side]).evalf(0, prec));
		}
		if ((fx[side]<0 && fx[!side]<0) || (fx[side]>0 && fx[!side]>0)) {
			// Oops, the root isn't bracketed any more.
			// Restore, and perform a bisection!
			xx[side] = xxprev;
			fx[side] = fxprev;

			// Ah, the bisection! Bisections converge linearly. Unfortunately,
			// they occur pretty often when Newton-Raphson arrives at an x too
			// close to the result on one side of the interval and
			// f(x-f(x)/f'(x)) turns out to have the same sign as f(x) due to
			// precision errors! Recall that this function does not have a
			// precision goal as one of its arguments but instead relies on
			// x converging to a fixed point. We speed up the (safe but slow)
			// bisection method by mixing in a dash of the (unsafer but faster)
			// secant method: Instead of splitting the interval at the
			// arithmetic mean (bisection), we split it nearer to the root as
			// determined by the secant between the values xx[0] and xx[1].
			// Don't set the secant_weight to one because that could disturb
			// the convergence in some corner cases!
			static const double secant_weight = 0.984375;  // == 63/64 < 1
			numeric xxmid = (1-secant_weight)*0.5*(xx[0]+xx[1])
			    + secant_weight*(xx[0]+fx[0]*(xx[0]-xx[1])/(fx[1]-fx[0]));
			numeric fxmid = ex_to<numeric>(f.subs(x==xxmid).evalf(0, prec));
			if (fxmid.is_zero()) {
				// Luck strikes...
				return xxmid;
			}
			if ((fxmid<0 && fx[side]>0) || (fxmid>0 && fx[side]<0)) {
				side = !side;
			}
			xxprev = xx[side];
			fxprev = fx[side];
			xx[side] = xxmid;
			fx[side] = fxmid;
		}
	} while (xxprev!=xx[side]);
	return xxprev;
}


/* Force inclusion of functions from inifcns_gamma and inifcns_zeta
 * for static lib (so ginsh will see them). */
unsigned force_include_tgamma = tgamma_SERIAL::serial;
unsigned force_include_zeta1 = zeta1_SERIAL::serial;

} // namespace GiNaC
