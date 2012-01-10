/** @file integral.cpp
 *
 *  Implementation of GiNaC's symbolic  integral. */

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

#include "integral.h"
#include "numeric.h"
#include "symbol.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "inifcns.h"
#include "wildcard.h"
#include "archive.h"
#include "registrar.h"
#include "utils.h"
#include "operators.h"
#include "relational.h"

using namespace std;

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(integral, basic,
  print_func<print_dflt>(&integral::do_print).
  print_func<print_latex>(&integral::do_print_latex))


//////////
// default constructor
//////////

integral::integral()
		: inherited(&integral::tinfo_static),
		x((new symbol())->setflag(status_flags::dynallocated))
{}

//////////
// other constructors
//////////

// public

integral::integral(const ex & x_, const ex & a_, const ex & b_, const ex & f_)
		: inherited(&integral::tinfo_static), x(x_), a(a_), b(b_), f(f_)
{
	if (!is_a<symbol>(x)) {
		throw(std::invalid_argument("first argument of integral must be of type symbol"));
	}
}

//////////
// archiving
//////////

integral::integral(const archive_node & n, lst & sym_lst) : inherited(n, sym_lst)
{
	n.find_ex("x", x, sym_lst);
	n.find_ex("a", a, sym_lst);
	n.find_ex("b", b, sym_lst);
	n.find_ex("f", f, sym_lst);
}

void integral::archive(archive_node & n) const
{
	inherited::archive(n);
	n.add_ex("x", x);
	n.add_ex("a", a);
	n.add_ex("b", b);
	n.add_ex("f", f);
}

DEFAULT_UNARCHIVE(integral)

//////////
// functions overriding virtual functions from base classes
//////////

void integral::do_print(const print_context & c, unsigned level) const
{
	c.s << "integral(";
	x.print(c);
	c.s << ",";
	a.print(c);
	c.s << ",";
	b.print(c);
	c.s << ",";
	f.print(c);
	c.s << ")";
}

void integral::do_print_latex(const print_latex & c, unsigned level) const
{
	string varname = ex_to<symbol>(x).get_name();
	if (level > precedence())
		c.s << "\\left(";
	c.s << "\\int_{";
	a.print(c);
	c.s << "}^{";
	b.print(c);
	c.s << "} d";
	if (varname.size() > 1)
		c.s << "\\," << varname << "\\:";
	else
		c.s << varname << "\\,";
	f.print(c,precedence());
	if (level > precedence())
		c.s << "\\right)";
}

int integral::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<integral>(other));
	const integral &o = static_cast<const integral &>(other);

	int cmpval = x.compare(o.x);
	if (cmpval)
		return cmpval;
	cmpval = a.compare(o.a);
	if (cmpval)
		return cmpval;
	cmpval = b.compare(o.b);
	if (cmpval)
		return cmpval;
	return f.compare(o.f);
}

ex integral::eval(int level) const
{
	if ((level==1) && (flags & status_flags::evaluated))
		return *this;
	if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));

	ex eintvar = (level==1) ? x : x.eval(level-1);
	ex ea      = (level==1) ? a : a.eval(level-1);
	ex eb      = (level==1) ? b : b.eval(level-1);
	ex ef      = (level==1) ? f : f.eval(level-1);

	if (!ef.has(eintvar) && !haswild(ef))
		return eb*ef-ea*ef;

	if (ea==eb)
		return _ex0;

	if (are_ex_trivially_equal(eintvar,x) && are_ex_trivially_equal(ea,a)
			&& are_ex_trivially_equal(eb,b) && are_ex_trivially_equal(ef,f))
		return this->hold();
	return (new integral(eintvar, ea, eb, ef))
		->setflag(status_flags::dynallocated | status_flags::evaluated);
}

ex integral::evalf(int level) const
{
	ex ea;
	ex eb;
	ex ef;

	if (level==1) {
		ea = a;
		eb = b;
		ef = f;
	} else if (level == -max_recursion_level) {
		throw(runtime_error("max recursion level reached"));
	} else {
		ea = a.evalf(level-1);
		eb = b.evalf(level-1);
		ef = f.evalf(level-1);
	}

	// 12.34 is just an arbitrary number used to check whether a number
	// results after subsituting a number for the integration variable.
	if (is_exactly_a<numeric>(ea) && is_exactly_a<numeric>(eb) 
			&& is_exactly_a<numeric>(ef.subs(x==12.34).evalf())) {
		try {
			return adaptivesimpson(x, ea, eb, ef);
		} catch (runtime_error &rte) {}
	}

	if (are_ex_trivially_equal(a, ea) && are_ex_trivially_equal(b, eb)
				&& are_ex_trivially_equal(f, ef))
			return *this;
		else
			return (new integral(x, ea, eb, ef))
				->setflag(status_flags::dynallocated);
}

int integral::max_integration_level = 15;
ex integral::relative_integration_error = 1e-8;

ex subsvalue(const ex & var, const ex & value, const ex & fun)
{
	ex result = fun.subs(var==value).evalf();
	if (is_a<numeric>(result))
		return result;
	throw logic_error("integrand does not evaluate to numeric");
}

struct error_and_integral
{
	error_and_integral(const ex &err, const ex &integ)
		:error(err), integral(integ){}
	ex error;
	ex integral;
};

struct error_and_integral_is_less
{
	bool operator()(const error_and_integral &e1,const error_and_integral &e2) const
	{
		int c = e1.integral.compare(e2.integral);
		if(c < 0)
			return true;
		if(c > 0)
			return false;
		return ex_is_less()(e1.error, e2.error);
	}
};

typedef map<error_and_integral, ex, error_and_integral_is_less> lookup_map;

/** Numeric integration routine based upon the "Adaptive Quadrature" one
  * in "Numerical Analysis" by Burden and Faires. Parameters are integration
  * variable, left boundary, right boundary, function to be integrated and
  * the relative integration error. The function should evalf into a number
  * after substituting the integration variable by a number. Another thing
  * to note is that this implementation is no good at integrating functions
  * with discontinuities. */
ex adaptivesimpson(const ex & x, const ex & a_in, const ex & b_in, const ex & f, const ex & error)
{
	// Check whether boundaries and error are numbers.
	ex a = is_exactly_a<numeric>(a_in) ? a_in : a_in.evalf();
	ex b = is_exactly_a<numeric>(b_in) ? b_in : b_in.evalf();
	if(!is_exactly_a<numeric>(a) || !is_exactly_a<numeric>(b))
		throw std::runtime_error("For numerical integration the boundaries of the integral should evalf into numbers.");
	if(!is_exactly_a<numeric>(error))
		throw std::runtime_error("For numerical integration the error should be a number.");

	// Use lookup table to be potentially much faster.
	static lookup_map lookup;
	static symbol ivar("ivar");
	ex lookupex = integral(ivar,a,b,f.subs(x==ivar));
	lookup_map::iterator emi = lookup.find(error_and_integral(error, lookupex));
	if (emi!=lookup.end())
		return emi->second;

	ex app = 0;
	int i = 1;
	exvector avec(integral::max_integration_level+1);
	exvector hvec(integral::max_integration_level+1);
	exvector favec(integral::max_integration_level+1);
	exvector fbvec(integral::max_integration_level+1);
	exvector fcvec(integral::max_integration_level+1);
	exvector svec(integral::max_integration_level+1);
	exvector errorvec(integral::max_integration_level+1);
	vector<int> lvec(integral::max_integration_level+1);

	avec[i] = a;
	hvec[i] = (b-a)/2;
	favec[i] = subsvalue(x, a, f);
	fcvec[i] = subsvalue(x, a+hvec[i], f);
	fbvec[i] = subsvalue(x, b, f);
	svec[i] = hvec[i]*(favec[i]+4*fcvec[i]+fbvec[i])/3;
	lvec[i] = 1;
	errorvec[i] = error*abs(svec[i]);

	while (i>0) {
		ex fd = subsvalue(x, avec[i]+hvec[i]/2, f);
		ex fe = subsvalue(x, avec[i]+3*hvec[i]/2, f);
		ex s1 = hvec[i]*(favec[i]+4*fd+fcvec[i])/6;
		ex s2 = hvec[i]*(fcvec[i]+4*fe+fbvec[i])/6;
		ex nu1 = avec[i];
		ex nu2 = favec[i];
		ex nu3 = fcvec[i];
		ex nu4 = fbvec[i];
		ex nu5 = hvec[i];
		// hopefully prevents a crash if the function is zero sometimes.
		ex nu6 = max(errorvec[i], abs(s1+s2)*error);
		ex nu7 = svec[i];
		int nu8 = lvec[i];
		--i;
		if (abs(ex_to<numeric>(s1+s2-nu7)) <= nu6)
			app+=(s1+s2);
		else {
			if (nu8>=integral::max_integration_level)
				throw runtime_error("max integration level reached");
			++i;
			avec[i] = nu1+nu5;
			favec[i] = nu3;
			fcvec[i] = fe;
			fbvec[i] = nu4;
			hvec[i] = nu5/2;
			errorvec[i]=nu6/2;
			svec[i] = s2;
			lvec[i] = nu8+1;
			++i;
			avec[i] = nu1;
			favec[i] = nu2;
			fcvec[i] = fd;
			fbvec[i] = nu3;
			hvec[i] = hvec[i-1];
			errorvec[i]=errorvec[i-1];
			svec[i] = s1;
			lvec[i] = lvec[i-1];
		}
	}

	lookup[error_and_integral(error, lookupex)]=app;
	return app;
}

int integral::degree(const ex & s) const
{
	return ((b-a)*f).degree(s);
}

int integral::ldegree(const ex & s) const
{
	return ((b-a)*f).ldegree(s);
}

ex integral::eval_ncmul(const exvector & v) const
{
	return f.eval_ncmul(v);
}

size_t integral::nops() const
{
	return 4;
}

ex integral::op(size_t i) const
{
	GINAC_ASSERT(i<4);

	switch (i) {
		case 0:
			return x;
		case 1:
			return a;
		case 2:
			return b;
		case 3:
			return f;
		default:
			throw (std::out_of_range("integral::op() out of range"));
	}
}

ex & integral::let_op(size_t i)
{
	ensure_if_modifiable();
	switch (i) {
		case 0:
			return x;
		case 1:
			return a;
		case 2:
			return b;
		case 3:
			return f;
		default:
			throw (std::out_of_range("integral::let_op() out of range"));
	}
}

ex integral::expand(unsigned options) const
{
	if (options==0 && (flags & status_flags::expanded))
		return *this;

	ex newa = a.expand(options);
	ex newb = b.expand(options);
	ex newf = f.expand(options);

	if (is_a<add>(newf)) {
		exvector v;
		v.reserve(newf.nops());
		for (size_t i=0; i<newf.nops(); ++i)
			v.push_back(integral(x, newa, newb, newf.op(i)).expand(options));
		return ex(add(v)).expand(options);
	}

	if (is_a<mul>(newf)) {
		ex prefactor = 1;
		ex rest = 1;
		for (size_t i=0; i<newf.nops(); ++i)
			if (newf.op(i).has(x))
				rest *= newf.op(i);
			else
				prefactor *= newf.op(i);
		if (prefactor != 1)
			return (prefactor*integral(x, newa, newb, rest)).expand(options);
	}

	if (are_ex_trivially_equal(a, newa) && are_ex_trivially_equal(b, newb)
			&& are_ex_trivially_equal(f, newf)) {
		if (options==0)
			this->setflag(status_flags::expanded);
		return *this;
	}

	const basic & newint = (new integral(x, newa, newb, newf))
		->setflag(status_flags::dynallocated);
	if (options == 0)
		newint.setflag(status_flags::expanded);
	return newint;
}

ex integral::derivative(const symbol & s) const
{
	if (s==x)
		throw(logic_error("differentiation with respect to dummy variable"));
	return b.diff(s)*f.subs(x==b)-a.diff(s)*f.subs(x==a)+integral(x, a, b, f.diff(s));
}

unsigned integral::return_type() const
{
	return f.return_type();
}

tinfo_t integral::return_type_tinfo() const
{
	return f.return_type_tinfo();
}

ex integral::conjugate() const
{
	ex conja = a.conjugate();
	ex conjb = b.conjugate();
	ex conjf = f.conjugate().subs(x.conjugate()==x);

	if (are_ex_trivially_equal(a, conja) && are_ex_trivially_equal(b, conjb)
			&& are_ex_trivially_equal(f, conjf))
		return *this;

	return (new integral(x, conja, conjb, conjf))
		->setflag(status_flags::dynallocated);
}

ex integral::eval_integ() const
{
	if (!(flags & status_flags::expanded))
		return this->expand().eval_integ();
	
	if (f==x)
		return b*b/2-a*a/2;
	if (is_a<power>(f) && f.op(0)==x) {
		if (f.op(1)==-1)
			return log(b/a);
		if (!f.op(1).has(x)) {
			ex primit = power(x,f.op(1)+1)/(f.op(1)+1);
			return primit.subs(x==b)-primit.subs(x==a);
		}
	}

	return *this;
}

} // namespace GiNaC
