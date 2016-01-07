#include "inifcns.h"

#include "infinity.h"
#include "lst.h"
#include "constant.h"
#include "numeric.h"
#include "operators.h"
#include "power.h"
#include "pseries.h"
#include "utils.h"

#include <sstream>
#include <stdexcept>
#include <vector>


namespace GiNaC {

// Stieltjes constants

static ex stieltjes1_evalf(const ex& x, PyObject* parent)
{
        if (is_exactly_a<numeric>(x)) {
		try {
			return stieltjes(ex_to<numeric>(x.evalf(0, parent)));
		} catch (const dunno &e) { }
	}

	return stieltjes(x).hold();
}


static ex stieltjes1_eval(const ex& arg)
{
	if (not arg.info(info_flags::numeric))
                return stieltjes(arg).hold();

	if (not arg.info(info_flags::integer))
                return stieltjes(ex_to<numeric>(arg));

        if (ex_to<numeric>(arg).is_zero())
                return Euler;
        else if (not ex_to<numeric>(arg).is_negative())
                return stieltjes(arg).hold();
        else
                throw std::runtime_error("Stieltjes constant of negative index");
}

static void stieltjes1_print_latex(const ex& arg, const print_context& c)
{
	c.s << "\\gamma";
        if (is_exactly_a<numeric>(arg) and ex_to<numeric>(arg).is_zero())
                return;
        c.s << "_{";
	arg.print(c);
	c.s << "}";
}


unsigned stieltjes1_SERIAL::serial = function::register_new(function_options("stieltjes", 1).
                                evalf_func(stieltjes1_evalf).
                                eval_func(stieltjes1_eval).
                                print_func<print_latex>(stieltjes1_print_latex).
                                overloaded(2));

//////////////////////////////////////////////////////////////////////
//
// Multiple zeta values  zeta(x)
//
// GiNaC function
//
//////////////////////////////////////////////////////////////////////


static ex zeta1_evalf(const ex& x, PyObject* parent)
{
        /*
	if (is_exactly_a<lst>(x) && (x.nops()>1)) {

		// multiple zeta value
		const int count = x.nops();
		const lst& xlst = ex_to<lst>(x);
		std::vector<int> r(count);

		// check parameters and convert them
		auto it1 = xlst.begin();
		auto it2 = r.begin();
		do {
			if (!(*it1).info(info_flags::posint)) {
				return zeta(x).hold();
			}
			*it2 = ex_to<numeric>(*it1).to_int();
			it1++;
			it2++;
		} while (it2 != r.end());

		// check for divergence
		if (r[0] == 1) {
			return zeta(x).hold();
		}

		// decide on summation algorithm
		// this is still a bit clumsy
		int limit = 10;
		if ((r[0] < limit) || ((count > 3) && (r[1] < limit/2))) {
			return numeric(zeta_do_sum_Crandall(r));
		} else {
			return numeric(zeta_do_sum_simple(r));
		}
	}*/

	// single zeta value
	if (x == 1) {
		return UnsignedInfinity;
	} else	if (is_exactly_a<numeric>(x)) {
		try {
			return zeta(ex_to<numeric>(x.evalf(0, parent)));
		} catch (const dunno &e) { }
	}

	return zeta(x).hold();
}


static ex zeta1_eval(const ex& m)
{
	if (is_exactly_a<lst>(m)) {
		if (m.nops() == 1) {
			return zeta(m.op(0));
		}
		return zeta(m).hold();
	}

	if (m.info(info_flags::numeric)) {
		const numeric& y = ex_to<numeric>(m);
		// trap integer arguments:
		if (y.is_integer()) {
			if (y.is_zero()) {
				return _ex_1_2;
			}
			if (y.is_equal(*_num1_p)) {
				//return zeta(m).hold();
				return UnsignedInfinity;
			}
			if (y.info(info_flags::posint)) {
				if (y.info(info_flags::odd)) {
					return zeta(m).hold();
				} else {
					return abs(bernoulli(y)) * pow(Pi, y) * pow(*_num2_p, y-(*_num1_p)) / factorial(y);
				}
			} else {
				if (y.info(info_flags::odd)) {
					return -bernoulli((*_num1_p)-y) / ((*_num1_p)-y);
				} else {
					return _ex0;
				}
			}
		}
		// zeta(float)
		if (y.info(info_flags::numeric) && !y.info(info_flags::crational)) {
			return zeta(y); // y is numeric
		}
	}
	return zeta(m).hold();
}

static ex zeta1_series(const ex& m, const relational& rel, int order, unsigned options)
{
	// use taylor expansion everywhere except at the singularity at 1
	const numeric val = ex_to<numeric>(m.subs(rel, subs_options::no_pattern));
	if (val != 1)
		throw do_taylor();  // caught by function::series()
	// at 1, use zeta's functional equation and develop the resulting expression
	return (pow(2,m) * pow(Pi, m-1) * sin(Pi*m/2) * tgamma(1-m) * zeta(1-m)).series(rel, order, options);
}

static ex zeta1_deriv(const ex& m, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	if (is_exactly_a<lst>(m)) {
		return _ex0;
	} else {
		return zetaderiv(_ex1, m);
	}
}


static void zeta1_print_latex(const ex& m_, const print_context& c)
{
	c.s << "\\zeta(";
	if (is_a<lst>(m_)) {
		const lst& m = ex_to<lst>(m_);
		auto it = m.begin();
		(*it).print(c);
		it++;
		for (; it != m.end(); it++) {
			c.s << ",";
			(*it).print(c);
		}
	} else {
		m_.print(c);
	}
	c.s << ")";
}


unsigned zeta1_SERIAL::serial = function::register_new(function_options("zeta", 1).
                                evalf_func(zeta1_evalf).
                                eval_func(zeta1_eval).
                                derivative_func(zeta1_deriv).
                                series_func(zeta1_series).
                                print_func<print_latex>(zeta1_print_latex).
                                overloaded(2));


//////////////////////////////////////////////////////////////////////
//
// Alternating Euler sum  zeta(x,s)
//
// GiNaC function
//
//////////////////////////////////////////////////////////////////////


static ex zeta2_evalf(const ex& x, const ex& s, PyObject* parent)
{
        /*
	if (is_exactly_a<lst>(x)) {

		// alternating Euler sum
		const int count = x.nops();
		const lst& xlst = ex_to<lst>(x);
		const lst& slst = ex_to<lst>(s);
		std::vector<int> xi(count);
		std::vector<int> si(count);

		// check parameters and convert them
		auto it_xread = xlst.begin();
		auto it_sread = slst.begin();
		auto it_xwrite = xi.begin();
		auto it_swrite = si.begin();
		do {
			if (!(*it_xread).info(info_flags::posint)) {
				return zeta(x, s).hold();
			}
			*it_xwrite = ex_to<numeric>(*it_xread).to_int();
			if (*it_sread > 0) {
				*it_swrite = 1;
			} else {
				*it_swrite = -1;
			}
			it_xread++;
			it_sread++;
			it_xwrite++;
			it_swrite++;
		} while (it_xwrite != xi.end());

		// check for divergence
		if ((xi[0] == 1) && (si[0] == 1)) {
			return zeta(x, s).hold();
		}

		// use Hoelder convolution
		return numeric(zeta_do_Hoelder_convolution(xi, si));
	}*/

	return zeta(x, s).hold();
}


static ex zeta2_eval(const ex& m, const ex& s_)
{
	if (is_exactly_a<lst>(s_)) {
		const lst& s = ex_to<lst>(s_);
		for (const auto & elem : s) {
			if ((elem).info(info_flags::positive)) {
				continue;
			}
			return zeta(m, s_).hold();
		}
		return zeta(m);
	} else if (s_.info(info_flags::positive)) {
		return zeta(m);
	}

	return zeta(m, s_).hold();
}


static ex zeta2_deriv(const ex& m, const ex& s, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	if (is_exactly_a<lst>(m)) {
		return _ex0;
	} else {
		if ((is_exactly_a<lst>(s) && s.op(0).info(info_flags::positive)) || s.info(info_flags::positive)) {
			return zetaderiv(_ex1, m);
		}
		return _ex0;
	}
}


static void zeta2_print_latex(const ex& m_, const ex& s_, const print_context& c)
{
	lst m;
	if (is_a<lst>(m_)) {
		m = ex_to<lst>(m_);
	} else {
		m = lst(m_);
	}
	lst s;
	if (is_a<lst>(s_)) {
		s = ex_to<lst>(s_);
	} else {
		s = lst(s_);
	}
	c.s << "\\zeta(";
	auto itm = m.begin();
	auto its = s.begin();
	if (*its < 0) {
		c.s << "\\overline{";
		(*itm).print(c);
		c.s << "}";
	} else {
		(*itm).print(c);
	}
	its++;
	itm++;
	for (; itm != m.end(); itm++, its++) {
		c.s << ",";
		if (*its < 0) {
			c.s << "\\overline{";
			(*itm).print(c);
			c.s << "}";
		} else {
			(*itm).print(c);
		}
	}
	c.s << ")";
}


unsigned zeta2_SERIAL::serial = function::register_new(function_options("zeta", 2).
                                evalf_func(zeta2_evalf).
                                eval_func(zeta2_eval).
                                derivative_func(zeta2_deriv).
                                print_func<print_latex>(zeta2_print_latex).
                                overloaded(2));

} // namespace GiNaC
