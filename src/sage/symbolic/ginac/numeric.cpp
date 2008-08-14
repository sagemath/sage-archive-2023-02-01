/** @file numeric.cpp
 *
 *  This file contains the interface to the underlying bignum package.
 *  Its most important design principle is to completely hide the inner
 *  working of that other package from the user of GiNaC.  It must either 
 *  provide implementation of arithmetic operators and numerical evaluation
 *  of special functions or implement the interface to the bignum package. */

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

#include "config.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <sstream>
#include <limits>

#include "numeric.h"
#include "ex.h"
#include "operators.h"
#include "archive.h"
#include "tostring.h"
#include "utils.h"

// CLN should pollute the global namespace as little as possible.  Hence, we
// include most of it here and include only the part needed for properly
// declaring cln::cl_number in numeric.h.  This can only be safely done in
// namespaced versions of CLN, i.e. version > 1.1.0.  Also, we only need a
// subset of CLN, so we don't include the complete <cln/cln.h> but only the
// essential stuff:
#include <cln/output.h>
#include <cln/integer_io.h>
#include <cln/integer_ring.h>
#include <cln/rational_io.h>
#include <cln/rational_ring.h>
#include <cln/lfloat_class.h>
#include <cln/lfloat_io.h>
#include <cln/real_io.h>
#include <cln/real_ring.h>
#include <cln/complex_io.h>
#include <cln/complex_ring.h>
#include <cln/numtheory.h>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(numeric, basic,
  print_func<print_context>(&numeric::do_print).
  print_func<print_latex>(&numeric::do_print_latex).
  print_func<print_csrc>(&numeric::do_print_csrc).
  print_func<print_csrc_cl_N>(&numeric::do_print_csrc_cl_N).
  print_func<print_tree>(&numeric::do_print_tree).
  print_func<print_python_repr>(&numeric::do_print_python_repr))

//////////
// default constructor
//////////

/** default ctor. Numerically it initializes to an integer zero. */
numeric::numeric() : basic(&numeric::tinfo_static)
{
	value = cln::cl_I(0);
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

// public

numeric::numeric(int i) : basic(&numeric::tinfo_static)
{
	// Not the whole int-range is available if we don't cast to long
	// first.  This is due to the behaviour of the cl_I-ctor, which
	// emphasizes efficiency.  However, if the integer is small enough
	// we save space and dereferences by using an immediate type.
	// (C.f. <cln/object.h>)
	// The #if clause prevents compiler warnings on 64bit machines where the
	// comparision is always true.
#if cl_value_len >= 32
	value = cln::cl_I(i);
#else
	if (i < (1L << (cl_value_len-1)) && i >= -(1L << (cl_value_len-1)))
		value = cln::cl_I(i);
	else
		value = cln::cl_I(static_cast<long>(i));
#endif
	setflag(status_flags::evaluated | status_flags::expanded);
}


numeric::numeric(unsigned int i) : basic(&numeric::tinfo_static)
{
	// Not the whole uint-range is available if we don't cast to ulong
	// first.  This is due to the behaviour of the cl_I-ctor, which
	// emphasizes efficiency.  However, if the integer is small enough
	// we save space and dereferences by using an immediate type.
	// (C.f. <cln/object.h>)
	// The #if clause prevents compiler warnings on 64bit machines where the
	// comparision is always true.
#if cl_value_len >= 32
	value = cln::cl_I(i);
#else
	if (i < (1UL << (cl_value_len-1)))
		value = cln::cl_I(i);
	else
		value = cln::cl_I(static_cast<unsigned long>(i));
#endif
	setflag(status_flags::evaluated | status_flags::expanded);
}


numeric::numeric(long i) : basic(&numeric::tinfo_static)
{
	value = cln::cl_I(i);
	setflag(status_flags::evaluated | status_flags::expanded);
}


numeric::numeric(unsigned long i) : basic(&numeric::tinfo_static)
{
	value = cln::cl_I(i);
	setflag(status_flags::evaluated | status_flags::expanded);
}


/** Constructor for rational numerics a/b.
 *
 *  @exception overflow_error (division by zero) */
numeric::numeric(long numer, long denom) : basic(&numeric::tinfo_static)
{
	if (!denom)
		throw std::overflow_error("division by zero");
	value = cln::cl_I(numer) / cln::cl_I(denom);
	setflag(status_flags::evaluated | status_flags::expanded);
}


numeric::numeric(double d) : basic(&numeric::tinfo_static)
{
	// We really want to explicitly use the type cl_LF instead of the
	// more general cl_F, since that would give us a cl_DF only which
	// will not be promoted to cl_LF if overflow occurs:
	value = cln::cl_float(d, cln::default_float_format);
	setflag(status_flags::evaluated | status_flags::expanded);
}


/** ctor from C-style string.  It also accepts complex numbers in GiNaC
 *  notation like "2+5*I". */
numeric::numeric(const char *s) : basic(&numeric::tinfo_static)
{
	cln::cl_N ctorval = 0;
	// parse complex numbers (functional but not completely safe, unfortunately
	// std::string does not understand regexpese):
	// ss should represent a simple sum like 2+5*I
	std::string ss = s;
	std::string::size_type delim;

	// make this implementation safe by adding explicit sign
	if (ss.at(0) != '+' && ss.at(0) != '-' && ss.at(0) != '#')
		ss = '+' + ss;

	// We use 'E' as exponent marker in the output, but some people insist on
	// writing 'e' at input, so let's substitute them right at the beginning:
	while ((delim = ss.find("e"))!=std::string::npos)
		ss.replace(delim,1,"E");

	// main parser loop:
	do {
		// chop ss into terms from left to right
		std::string term;
		bool imaginary = false;
		delim = ss.find_first_of(std::string("+-"),1);
		// Do we have an exponent marker like "31.415E-1"?  If so, hop on!
		if (delim!=std::string::npos && ss.at(delim-1)=='E')
			delim = ss.find_first_of(std::string("+-"),delim+1);
		term = ss.substr(0,delim);
		if (delim!=std::string::npos)
			ss = ss.substr(delim);
		// is the term imaginary?
		if (term.find("I")!=std::string::npos) {
			// erase 'I':
			term.erase(term.find("I"),1);
			// erase '*':
			if (term.find("*")!=std::string::npos)
				term.erase(term.find("*"),1);
			// correct for trivial +/-I without explicit factor on I:
			if (term.size()==1)
				term += '1';
			imaginary = true;
		}
		if (term.find('.')!=std::string::npos || term.find('E')!=std::string::npos) {
			// CLN's short type cl_SF is not very useful within the GiNaC
			// framework where we are mainly interested in the arbitrary
			// precision type cl_LF.  Hence we go straight to the construction
			// of generic floats.  In order to create them we have to convert
			// our own floating point notation used for output and construction
			// from char * to CLN's generic notation:
			// 3.14      -->   3.14e0_<Digits>
			// 31.4E-1   -->   31.4e-1_<Digits>
			// and s on.
			// No exponent marker?  Let's add a trivial one.
			if (term.find("E")==std::string::npos)
				term += "E0";
			// E to lower case
			term = term.replace(term.find("E"),1,"e");
			// append _<Digits> to term
			term += "_" + ToString((unsigned)Digits);
			// construct float using cln::cl_F(const char *) ctor.
			if (imaginary)
				ctorval = ctorval + cln::complex(cln::cl_I(0),cln::cl_F(term.c_str()));
			else
				ctorval = ctorval + cln::cl_F(term.c_str());
		} else {
			// this is not a floating point number...
			if (imaginary)
				ctorval = ctorval + cln::complex(cln::cl_I(0),cln::cl_R(term.c_str()));
			else
				ctorval = ctorval + cln::cl_R(term.c_str());
		}
	} while (delim != std::string::npos);
	value = ctorval;
	setflag(status_flags::evaluated | status_flags::expanded);
}


/** Ctor from CLN types.  This is for the initiated user or internal use
 *  only. */
numeric::numeric(const cln::cl_N &z) : basic(&numeric::tinfo_static)
{
	value = z;
	setflag(status_flags::evaluated | status_flags::expanded);
}


//////////
// archiving
//////////

numeric::numeric(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	cln::cl_N ctorval = 0;

	// Read number as string
	std::string str;
	if (n.find_string("number", str)) {
		std::istringstream s(str);
		cln::cl_idecoded_float re, im;
		char c;
		s.get(c);
		switch (c) {
			case 'R':    // Integer-decoded real number
				s >> re.sign >> re.mantissa >> re.exponent;
				ctorval = re.sign * re.mantissa * cln::expt(cln::cl_float(2.0, cln::default_float_format), re.exponent);
				break;
			case 'C':    // Integer-decoded complex number
				s >> re.sign >> re.mantissa >> re.exponent;
				s >> im.sign >> im.mantissa >> im.exponent;
				ctorval = cln::complex(re.sign * re.mantissa * cln::expt(cln::cl_float(2.0, cln::default_float_format), re.exponent),
				                       im.sign * im.mantissa * cln::expt(cln::cl_float(2.0, cln::default_float_format), im.exponent));
				break;
			default:    // Ordinary number
				s.putback(c);
				s >> ctorval;
				break;
		}
	}
	value = ctorval;
	setflag(status_flags::evaluated | status_flags::expanded);
}

void numeric::archive(archive_node &n) const
{
	inherited::archive(n);

	// Write number as string
	std::ostringstream s;
	if (this->is_crational())
		s << value;
	else {
		// Non-rational numbers are written in an integer-decoded format
		// to preserve the precision
		if (this->is_real()) {
			cln::cl_idecoded_float re = cln::integer_decode_float(cln::the<cln::cl_F>(value));
			s << "R";
			s << re.sign << " " << re.mantissa << " " << re.exponent;
		} else {
			cln::cl_idecoded_float re = cln::integer_decode_float(cln::the<cln::cl_F>(cln::realpart(cln::the<cln::cl_N>(value))));
			cln::cl_idecoded_float im = cln::integer_decode_float(cln::the<cln::cl_F>(cln::imagpart(cln::the<cln::cl_N>(value))));
			s << "C";
			s << re.sign << " " << re.mantissa << " " << re.exponent << " ";
			s << im.sign << " " << im.mantissa << " " << im.exponent;
		}
	}
	n.add_string("number", s.str());
}

DEFAULT_UNARCHIVE(numeric)

//////////
// functions overriding virtual functions from base classes
//////////

/** Helper function to print a real number in a nicer way than is CLN's
 *  default.  Instead of printing 42.0L0 this just prints 42.0 to ostream os
 *  and instead of 3.99168L7 it prints 3.99168E7.  This is fine in GiNaC as
 *  long as it only uses cl_LF and no other floating point types that we might
 *  want to visibly distinguish from cl_LF.
 *
 *  @see numeric::print() */
static void print_real_number(const print_context & c, const cln::cl_R & x)
{
	cln::cl_print_flags ourflags;
	if (cln::instanceof(x, cln::cl_RA_ring)) {
		// case 1: integer or rational
		if (cln::instanceof(x, cln::cl_I_ring) ||
		    !is_a<print_latex>(c)) {
			cln::print_real(c.s, ourflags, x);
		} else {  // rational output in LaTeX context
			if (x < 0)
				c.s << "-";
			c.s << "\\frac{";
			cln::print_real(c.s, ourflags, cln::abs(cln::numerator(cln::the<cln::cl_RA>(x))));
			c.s << "}{";
			cln::print_real(c.s, ourflags, cln::denominator(cln::the<cln::cl_RA>(x)));
			c.s << '}';
		}
	} else {
		// case 2: float
		// make CLN believe this number has default_float_format, so it prints
		// 'E' as exponent marker instead of 'L':
		ourflags.default_float_format = cln::float_format(cln::the<cln::cl_F>(x));
		cln::print_real(c.s, ourflags, x);
	}
}

/** Helper function to print integer number in C++ source format.
 *
 *  @see numeric::print() */
static void print_integer_csrc(const print_context & c, const cln::cl_I & x)
{
	// Print small numbers in compact float format, but larger numbers in
	// scientific format
	const int max_cln_int = 536870911; // 2^29-1
	if (x >= cln::cl_I(-max_cln_int) && x <= cln::cl_I(max_cln_int))
		c.s << cln::cl_I_to_int(x) << ".0";
	else
		c.s << cln::double_approx(x);
}

/** Helper function to print real number in C++ source format.
 *
 *  @see numeric::print() */
static void print_real_csrc(const print_context & c, const cln::cl_R & x)
{
	if (cln::instanceof(x, cln::cl_I_ring)) {

		// Integer number
		print_integer_csrc(c, cln::the<cln::cl_I>(x));

	} else if (cln::instanceof(x, cln::cl_RA_ring)) {

		// Rational number
		const cln::cl_I numer = cln::numerator(cln::the<cln::cl_RA>(x));
		const cln::cl_I denom = cln::denominator(cln::the<cln::cl_RA>(x));
		if (cln::plusp(x) > 0) {
			c.s << "(";
			print_integer_csrc(c, numer);
		} else {
			c.s << "-(";
			print_integer_csrc(c, -numer);
		}
		c.s << "/";
		print_integer_csrc(c, denom);
		c.s << ")";

	} else {

		// Anything else
		c.s << cln::double_approx(x);
	}
}

template<typename T1, typename T2> 
static inline bool coerce(T1& dst, const T2& arg);

/** 
 * @brief Check if CLN integer can be converted into int
 *
 * @sa http://www.ginac.de/pipermail/cln-list/2006-October/000248.html
 */
template<>
inline bool coerce<int, cln::cl_I>(int& dst, const cln::cl_I& arg)
{
	static const cln::cl_I cl_max_int = 
		(cln::cl_I)(long)(std::numeric_limits<int>::max());
	static const cln::cl_I cl_min_int =
		(cln::cl_I)(long)(std::numeric_limits<int>::min());
	if ((arg >= cl_min_int) && (arg <= cl_max_int)) {
		dst = cl_I_to_int(arg);
		return true;
	}
	return false;
}

template<>
inline bool coerce<unsigned int, cln::cl_I>(unsigned int& dst, const cln::cl_I& arg)
{
	static const cln::cl_I cl_max_uint = 
		(cln::cl_I)(unsigned long)(std::numeric_limits<unsigned int>::max());
	if ((! minusp(arg)) && (arg <= cl_max_uint)) {
		dst = cl_I_to_uint(arg);
		return true;
	}
	return false;
}

/** Helper function to print real number in C++ source format using cl_N types.
 *
 *  @see numeric::print() */
static void print_real_cl_N(const print_context & c, const cln::cl_R & x)
{
	if (cln::instanceof(x, cln::cl_I_ring)) {

		int dst;
		// fixnum
		if (coerce(dst, cln::the<cln::cl_I>(x))) {
			// can be converted to native int
			if (dst < 0)
				c.s << "(-" << dst << ")";
			else
				c.s << dst;
		} else {
			// bignum
			c.s << "cln::cl_I(\"";
			print_real_number(c, x);
			c.s << "\")";
		}
	} else if (cln::instanceof(x, cln::cl_RA_ring)) {

		// Rational number
		cln::cl_print_flags ourflags;
		c.s << "cln::cl_RA(\"";
		cln::print_rational(c.s, ourflags, cln::the<cln::cl_RA>(x));
		c.s << "\")";

	} else {

		// Anything else
		c.s << "cln::cl_F(\"";
		print_real_number(c, cln::cl_float(1.0, cln::default_float_format) * x);
		c.s << "_" << Digits << "\")";
	}
}

void numeric::print_numeric(const print_context & c, const char *par_open, const char *par_close, const char *imag_sym, const char *mul_sym, unsigned level) const
{
	const cln::cl_R r = cln::realpart(value);
	const cln::cl_R i = cln::imagpart(value);

	if (cln::zerop(i)) {

		// case 1, real:  x  or  -x
		if ((precedence() <= level) && (!this->is_nonneg_integer())) {
			c.s << par_open;
			print_real_number(c, r);
			c.s << par_close;
		} else {
			print_real_number(c, r);
		}

	} else {
		if (cln::zerop(r)) {

			// case 2, imaginary:  y*I  or  -y*I
			if (i == 1)
				c.s << imag_sym;
			else {
				if (precedence()<=level)
					c.s << par_open;
				if (i == -1)
					c.s << "-" << imag_sym;
				else {
					print_real_number(c, i);
					c.s << mul_sym << imag_sym;
				}
				if (precedence()<=level)
					c.s << par_close;
			}

		} else {

			// case 3, complex:  x+y*I  or  x-y*I  or  -x+y*I  or  -x-y*I
			if (precedence() <= level)
				c.s << par_open;
			print_real_number(c, r);
			if (i < 0) {
				if (i == -1) {
					c.s << "-" << imag_sym;
				} else {
					print_real_number(c, i);
					c.s << mul_sym << imag_sym;
				}
			} else {
				if (i == 1) {
					c.s << "+" << imag_sym;
				} else {
					c.s << "+";
					print_real_number(c, i);
					c.s << mul_sym << imag_sym;
				}
			}
			if (precedence() <= level)
				c.s << par_close;
		}
	}
}

void numeric::do_print(const print_context & c, unsigned level) const
{
	print_numeric(c, "(", ")", "I", "*", level);
}

void numeric::do_print_latex(const print_latex & c, unsigned level) const
{
	print_numeric(c, "{(", ")}", "i", " ", level);
}

void numeric::do_print_csrc(const print_csrc & c, unsigned level) const
{
	std::ios::fmtflags oldflags = c.s.flags();
	c.s.setf(std::ios::scientific);
	int oldprec = c.s.precision();

	// Set precision
	if (is_a<print_csrc_double>(c))
		c.s.precision(std::numeric_limits<double>::digits10 + 1);
	else
		c.s.precision(std::numeric_limits<float>::digits10 + 1);

	if (this->is_real()) {

		// Real number
		print_real_csrc(c, cln::the<cln::cl_R>(value));

	} else {

		// Complex number
		c.s << "std::complex<";
		if (is_a<print_csrc_double>(c))
			c.s << "double>(";
		else
			c.s << "float>(";

		print_real_csrc(c, cln::realpart(value));
		c.s << ",";
		print_real_csrc(c, cln::imagpart(value));
		c.s << ")";
	}

	c.s.flags(oldflags);
	c.s.precision(oldprec);
}

void numeric::do_print_csrc_cl_N(const print_csrc_cl_N & c, unsigned level) const
{
	if (this->is_real()) {

		// Real number
		print_real_cl_N(c, cln::the<cln::cl_R>(value));

	} else {

		// Complex number
		c.s << "cln::complex(";
		print_real_cl_N(c, cln::realpart(value));
		c.s << ",";
		print_real_cl_N(c, cln::imagpart(value));
		c.s << ")";
	}
}

void numeric::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << value
	    << " (" << class_name() << ")" << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << std::endl;
}

void numeric::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << "('";
	print_numeric(c, "(", ")", "I", "*", level);
	c.s << "')";
}

bool numeric::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::numeric:
		case info_flags::polynomial:
		case info_flags::rational_function:
		case info_flags::expanded:
			return true;
		case info_flags::real:
			return is_real();
		case info_flags::rational:
		case info_flags::rational_polynomial:
			return is_rational();
		case info_flags::crational:
		case info_flags::crational_polynomial:
			return is_crational();
		case info_flags::integer:
		case info_flags::integer_polynomial:
			return is_integer();
		case info_flags::cinteger:
		case info_flags::cinteger_polynomial:
			return is_cinteger();
		case info_flags::positive:
			return is_positive();
		case info_flags::negative:
			return is_negative();
		case info_flags::nonnegative:
			return !is_negative();
		case info_flags::posint:
			return is_pos_integer();
		case info_flags::negint:
			return is_integer() && is_negative();
		case info_flags::nonnegint:
			return is_nonneg_integer();
		case info_flags::even:
			return is_even();
		case info_flags::odd:
			return is_odd();
		case info_flags::prime:
			return is_prime();
		case info_flags::algebraic:
			return !is_real();
	}
	return false;
}

bool numeric::is_polynomial(const ex & var) const
{
	return true;
}

int numeric::degree(const ex & s) const
{
	return 0;
}

int numeric::ldegree(const ex & s) const
{
	return 0;
}

ex numeric::coeff(const ex & s, int n) const
{
	return n==0 ? *this : _ex0;
}

/** Disassemble real part and imaginary part to scan for the occurrence of a
 *  single number.  Also handles the imaginary unit.  It ignores the sign on
 *  both this and the argument, which may lead to what might appear as funny
 *  results:  (2+I).has(-2) -> true.  But this is consistent, since we also
 *  would like to have (-2+I).has(2) -> true and we want to think about the
 *  sign as a multiplicative factor. */
bool numeric::has(const ex &other, unsigned options) const
{
	if (!is_exactly_a<numeric>(other))
		return false;
	const numeric &o = ex_to<numeric>(other);
	if (this->is_equal(o) || this->is_equal(-o))
		return true;
	if (o.imag().is_zero()) {   // e.g. scan for 3 in -3*I
		if (!this->real().is_equal(*_num0_p))
			if (this->real().is_equal(o) || this->real().is_equal(-o))
				return true;
		if (!this->imag().is_equal(*_num0_p))
			if (this->imag().is_equal(o) || this->imag().is_equal(-o))
				return true;
		return false;
	}
	else {
		if (o.is_equal(I))  // e.g scan for I in 42*I
			return !this->is_real();
		if (o.real().is_zero())  // e.g. scan for 2*I in 2*I+1
			if (!this->imag().is_equal(*_num0_p))
				if (this->imag().is_equal(o*I) || this->imag().is_equal(-o*I))
					return true;
	}
	return false;
}


/** Evaluation of numbers doesn't do anything at all. */
ex numeric::eval(int level) const
{
	// Warning: if this is ever gonna do something, the ex ctors from all kinds
	// of numbers should be checking for status_flags::evaluated.
	return this->hold();
}


/** Cast numeric into a floating-point object.  For example exact numeric(1) is
 *  returned as a 1.0000000000000000000000 and so on according to how Digits is
 *  currently set.  In case the object already was a floating point number the
 *  precision is trimmed to match the currently set default.
 *
 *  @param level  ignored, only needed for overriding basic::evalf.
 *  @return  an ex-handle to a numeric. */
ex numeric::evalf(int level) const
{
	// level can safely be discarded for numeric objects.
	return numeric(cln::cl_float(1.0, cln::default_float_format) * value);
}

ex numeric::conjugate() const
{
	if (is_real()) {
		return *this;
	}
	return numeric(cln::conjugate(this->value));
}

ex numeric::real_part() const
{
	return numeric(cln::realpart(value));
}

ex numeric::imag_part() const
{
	return numeric(cln::imagpart(value));
}

// protected

int numeric::compare_same_type(const basic &other) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(other));
	const numeric &o = static_cast<const numeric &>(other);
	
	return this->compare(o);
}


bool numeric::is_equal_same_type(const basic &other) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(other));
	const numeric &o = static_cast<const numeric &>(other);
	
	return this->is_equal(o);
}


unsigned numeric::calchash() const
{
	// Base computation of hashvalue on CLN's hashcode.  Note: That depends
	// only on the number's value, not its type or precision (i.e. a true
	// equivalence relation on numbers).  As a consequence, 3 and 3.0 share
	// the same hashvalue.  That shouldn't really matter, though.
	setflag(status_flags::hash_calculated);
	hashvalue = golden_ratio_hash(cln::equal_hashcode(value));
	return hashvalue;
}


//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// public

/** Numerical addition method.  Adds argument to *this and returns result as
 *  a numeric object. */
const numeric numeric::add(const numeric &other) const
{
	return numeric(value + other.value);
}


/** Numerical subtraction method.  Subtracts argument from *this and returns
 *  result as a numeric object. */
const numeric numeric::sub(const numeric &other) const
{
	return numeric(value - other.value);
}


/** Numerical multiplication method.  Multiplies *this and argument and returns
 *  result as a numeric object. */
const numeric numeric::mul(const numeric &other) const
{
	return numeric(value * other.value);
}


/** Numerical division method.  Divides *this by argument and returns result as
 *  a numeric object.
 *
 *  @exception overflow_error (division by zero) */
const numeric numeric::div(const numeric &other) const
{
	if (cln::zerop(other.value))
		throw std::overflow_error("numeric::div(): division by zero");
	return numeric(value / other.value);
}


/** Numerical exponentiation.  Raises *this to the power given as argument and
 *  returns result as a numeric object. */
const numeric numeric::power(const numeric &other) const
{
	// Shortcut for efficiency and numeric stability (as in 1.0 exponent):
	// trap the neutral exponent.
	if (&other==_num1_p || cln::equal(other.value,_num1_p->value))
		return *this;
	
	if (cln::zerop(value)) {
		if (cln::zerop(other.value))
			throw std::domain_error("numeric::eval(): pow(0,0) is undefined");
		else if (cln::zerop(cln::realpart(other.value)))
			throw std::domain_error("numeric::eval(): pow(0,I) is undefined");
		else if (cln::minusp(cln::realpart(other.value)))
			throw std::overflow_error("numeric::eval(): division by zero");
		else
			return *_num0_p;
	}
	return numeric(cln::expt(value, other.value));
}



/** Numerical addition method.  Adds argument to *this and returns result as
 *  a numeric object on the heap.  Use internally only for direct wrapping into
 *  an ex object, where the result would end up on the heap anyways. */
const numeric &numeric::add_dyn(const numeric &other) const
{
	// Efficiency shortcut: trap the neutral element by pointer.  This hack
	// is supposed to keep the number of distinct numeric objects low.
	if (this==_num0_p)
		return other;
	else if (&other==_num0_p)
		return *this;
	
	return static_cast<const numeric &>((new numeric(value + other.value))->
	                                    setflag(status_flags::dynallocated));
}


/** Numerical subtraction method.  Subtracts argument from *this and returns
 *  result as a numeric object on the heap.  Use internally only for direct
 *  wrapping into an ex object, where the result would end up on the heap
 *  anyways. */
const numeric &numeric::sub_dyn(const numeric &other) const
{
	// Efficiency shortcut: trap the neutral exponent (first by pointer).  This
	// hack is supposed to keep the number of distinct numeric objects low.
	if (&other==_num0_p || cln::zerop(other.value))
		return *this;
	
	return static_cast<const numeric &>((new numeric(value - other.value))->
	                                    setflag(status_flags::dynallocated));
}


/** Numerical multiplication method.  Multiplies *this and argument and returns
 *  result as a numeric object on the heap.  Use internally only for direct
 *  wrapping into an ex object, where the result would end up on the heap
 *  anyways. */
const numeric &numeric::mul_dyn(const numeric &other) const
{
	// Efficiency shortcut: trap the neutral element by pointer.  This hack
	// is supposed to keep the number of distinct numeric objects low.
	if (this==_num1_p)
		return other;
	else if (&other==_num1_p)
		return *this;
	
	return static_cast<const numeric &>((new numeric(value * other.value))->
	                                    setflag(status_flags::dynallocated));
}


/** Numerical division method.  Divides *this by argument and returns result as
 *  a numeric object on the heap.  Use internally only for direct wrapping
 *  into an ex object, where the result would end up on the heap
 *  anyways.
 *
 *  @exception overflow_error (division by zero) */
const numeric &numeric::div_dyn(const numeric &other) const
{
	// Efficiency shortcut: trap the neutral element by pointer.  This hack
	// is supposed to keep the number of distinct numeric objects low.
	if (&other==_num1_p)
		return *this;
	if (cln::zerop(cln::the<cln::cl_N>(other.value)))
		throw std::overflow_error("division by zero");
	return static_cast<const numeric &>((new numeric(value / other.value))->
	                                    setflag(status_flags::dynallocated));
}


/** Numerical exponentiation.  Raises *this to the power given as argument and
 *  returns result as a numeric object on the heap.  Use internally only for
 *  direct wrapping into an ex object, where the result would end up on the
 *  heap anyways. */
const numeric &numeric::power_dyn(const numeric &other) const
{
	// Efficiency shortcut: trap the neutral exponent (first try by pointer, then
	// try harder, since calls to cln::expt() below may return amazing results for
	// floating point exponent 1.0).
	if (&other==_num1_p || cln::equal(other.value, _num1_p->value))
		return *this;
	
	if (cln::zerop(value)) {
		if (cln::zerop(other.value))
			throw std::domain_error("numeric::eval(): pow(0,0) is undefined");
		else if (cln::zerop(cln::realpart(other.value)))
			throw std::domain_error("numeric::eval(): pow(0,I) is undefined");
		else if (cln::minusp(cln::realpart(other.value)))
			throw std::overflow_error("numeric::eval(): division by zero");
		else
			return *_num0_p;
	}
	return static_cast<const numeric &>((new numeric(cln::expt(value, other.value)))->
	                                     setflag(status_flags::dynallocated));
}


const numeric &numeric::operator=(int i)
{
	return operator=(numeric(i));
}


const numeric &numeric::operator=(unsigned int i)
{
	return operator=(numeric(i));
}


const numeric &numeric::operator=(long i)
{
	return operator=(numeric(i));
}


const numeric &numeric::operator=(unsigned long i)
{
	return operator=(numeric(i));
}


const numeric &numeric::operator=(double d)
{
	return operator=(numeric(d));
}


const numeric &numeric::operator=(const char * s)
{
	return operator=(numeric(s));
}


/** Inverse of a number. */
const numeric numeric::inverse() const
{
	if (cln::zerop(value))
		throw std::overflow_error("numeric::inverse(): division by zero");
	return numeric(cln::recip(value));
}

/** Return the step function of a numeric. The imaginary part of it is
 *  ignored because the step function is generally considered real but
 *  a numeric may develop a small imaginary part due to rounding errors.
 */
numeric numeric::step() const
{	cln::cl_R r = cln::realpart(value);
	if(cln::zerop(r))
		return numeric(1,2);
	if(cln::plusp(r))
		return 1;
	return 0;
}

/** Return the complex half-plane (left or right) in which the number lies.
 *  csgn(x)==0 for x==0, csgn(x)==1 for Re(x)>0 or Re(x)=0 and Im(x)>0,
 *  csgn(x)==-1 for Re(x)<0 or Re(x)=0 and Im(x)<0.
 *
 *  @see numeric::compare(const numeric &other) */
int numeric::csgn() const
{
	if (cln::zerop(value))
		return 0;
	cln::cl_R r = cln::realpart(value);
	if (!cln::zerop(r)) {
		if (cln::plusp(r))
			return 1;
		else
			return -1;
	} else {
		if (cln::plusp(cln::imagpart(value)))
			return 1;
		else
			return -1;
	}
}


/** This method establishes a canonical order on all numbers.  For complex
 *  numbers this is not possible in a mathematically consistent way but we need
 *  to establish some order and it ought to be fast.  So we simply define it
 *  to be compatible with our method csgn.
 *
 *  @return csgn(*this-other)
 *  @see numeric::csgn() */
int numeric::compare(const numeric &other) const
{
	// Comparing two real numbers?
	if (cln::instanceof(value, cln::cl_R_ring) &&
		cln::instanceof(other.value, cln::cl_R_ring))
		// Yes, so just cln::compare them
		return cln::compare(cln::the<cln::cl_R>(value), cln::the<cln::cl_R>(other.value));
	else {
		// No, first cln::compare real parts...
		cl_signean real_cmp = cln::compare(cln::realpart(value), cln::realpart(other.value));
		if (real_cmp)
			return real_cmp;
		// ...and then the imaginary parts.
		return cln::compare(cln::imagpart(value), cln::imagpart(other.value));
	}
}


bool numeric::is_equal(const numeric &other) const
{
	return cln::equal(value, other.value);
}


/** True if object is zero. */
bool numeric::is_zero() const
{
	return cln::zerop(value);
}


/** True if object is not complex and greater than zero. */
bool numeric::is_positive() const
{
	if (cln::instanceof(value, cln::cl_R_ring))  // real?
		return cln::plusp(cln::the<cln::cl_R>(value));
	return false;
}


/** True if object is not complex and less than zero. */
bool numeric::is_negative() const
{
	if (cln::instanceof(value, cln::cl_R_ring))  // real?
		return cln::minusp(cln::the<cln::cl_R>(value));
	return false;
}


/** True if object is a non-complex integer. */
bool numeric::is_integer() const
{
	return cln::instanceof(value, cln::cl_I_ring);
}


/** True if object is an exact integer greater than zero. */
bool numeric::is_pos_integer() const
{
	return (cln::instanceof(value, cln::cl_I_ring) && cln::plusp(cln::the<cln::cl_I>(value)));
}


/** True if object is an exact integer greater or equal zero. */
bool numeric::is_nonneg_integer() const
{
	return (cln::instanceof(value, cln::cl_I_ring) && !cln::minusp(cln::the<cln::cl_I>(value)));
}


/** True if object is an exact even integer. */
bool numeric::is_even() const
{
	return (cln::instanceof(value, cln::cl_I_ring) && cln::evenp(cln::the<cln::cl_I>(value)));
}


/** True if object is an exact odd integer. */
bool numeric::is_odd() const
{
	return (cln::instanceof(value, cln::cl_I_ring) && cln::oddp(cln::the<cln::cl_I>(value)));
}


/** Probabilistic primality test.
 *
 *  @return  true if object is exact integer and prime. */
bool numeric::is_prime() const
{
	return (cln::instanceof(value, cln::cl_I_ring)  // integer?
	     && cln::plusp(cln::the<cln::cl_I>(value))  // positive?
	     && cln::isprobprime(cln::the<cln::cl_I>(value)));
}


/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_rational() const
{
	return cln::instanceof(value, cln::cl_RA_ring);
}


/** True if object is a real integer, rational or float (but not complex). */
bool numeric::is_real() const
{
	return cln::instanceof(value, cln::cl_R_ring);
}


bool numeric::operator==(const numeric &other) const
{
	return cln::equal(value, other.value);
}


bool numeric::operator!=(const numeric &other) const
{
	return !cln::equal(value, other.value);
}


/** True if object is element of the domain of integers extended by I, i.e. is
 *  of the form a+b*I, where a and b are integers. */
bool numeric::is_cinteger() const
{
	if (cln::instanceof(value, cln::cl_I_ring))
		return true;
	else if (!this->is_real()) {  // complex case, handle n+m*I
		if (cln::instanceof(cln::realpart(value), cln::cl_I_ring) &&
		    cln::instanceof(cln::imagpart(value), cln::cl_I_ring))
			return true;
	}
	return false;
}


/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_crational() const
{
	if (cln::instanceof(value, cln::cl_RA_ring))
		return true;
	else if (!this->is_real()) {  // complex case, handle Q(i):
		if (cln::instanceof(cln::realpart(value), cln::cl_RA_ring) &&
		    cln::instanceof(cln::imagpart(value), cln::cl_RA_ring))
			return true;
	}
	return false;
}


/** Numerical comparison: less.
 *
 *  @exception invalid_argument (complex inequality) */ 
bool numeric::operator<(const numeric &other) const
{
	if (this->is_real() && other.is_real())
		return (cln::the<cln::cl_R>(value) < cln::the<cln::cl_R>(other.value));
	throw std::invalid_argument("numeric::operator<(): complex inequality");
}


/** Numerical comparison: less or equal.
 *
 *  @exception invalid_argument (complex inequality) */ 
bool numeric::operator<=(const numeric &other) const
{
	if (this->is_real() && other.is_real())
		return (cln::the<cln::cl_R>(value) <= cln::the<cln::cl_R>(other.value));
	throw std::invalid_argument("numeric::operator<=(): complex inequality");
}


/** Numerical comparison: greater.
 *
 *  @exception invalid_argument (complex inequality) */ 
bool numeric::operator>(const numeric &other) const
{
	if (this->is_real() && other.is_real())
		return (cln::the<cln::cl_R>(value) > cln::the<cln::cl_R>(other.value));
	throw std::invalid_argument("numeric::operator>(): complex inequality");
}


/** Numerical comparison: greater or equal.
 *
 *  @exception invalid_argument (complex inequality) */  
bool numeric::operator>=(const numeric &other) const
{
	if (this->is_real() && other.is_real())
		return (cln::the<cln::cl_R>(value) >= cln::the<cln::cl_R>(other.value));
	throw std::invalid_argument("numeric::operator>=(): complex inequality");
}


/** Converts numeric types to machine's int.  You should check with
 *  is_integer() if the number is really an integer before calling this method.
 *  You may also consider checking the range first. */
int numeric::to_int() const
{
	GINAC_ASSERT(this->is_integer());
	return cln::cl_I_to_int(cln::the<cln::cl_I>(value));
}


/** Converts numeric types to machine's long.  You should check with
 *  is_integer() if the number is really an integer before calling this method.
 *  You may also consider checking the range first. */
long numeric::to_long() const
{
	GINAC_ASSERT(this->is_integer());
	return cln::cl_I_to_long(cln::the<cln::cl_I>(value));
}


/** Converts numeric types to machine's double. You should check with is_real()
 *  if the number is really not complex before calling this method. */
double numeric::to_double() const
{
	GINAC_ASSERT(this->is_real());
	return cln::double_approx(cln::realpart(value));
}


/** Returns a new CLN object of type cl_N, representing the value of *this.
 *  This method may be used when mixing GiNaC and CLN in one project.
 */
cln::cl_N numeric::to_cl_N() const
{
	return value;
}


/** Real part of a number. */
const numeric numeric::real() const
{
	return numeric(cln::realpart(value));
}


/** Imaginary part of a number. */
const numeric numeric::imag() const
{
	return numeric(cln::imagpart(value));
}


/** Numerator.  Computes the numerator of rational numbers, rationalized
 *  numerator of complex if real and imaginary part are both rational numbers
 *  (i.e numer(4/3+5/6*I) == 8+5*I), the number carrying the sign in all other
 *  cases. */
const numeric numeric::numer() const
{
	if (cln::instanceof(value, cln::cl_I_ring))
		return numeric(*this);  // integer case
	
	else if (cln::instanceof(value, cln::cl_RA_ring))
		return numeric(cln::numerator(cln::the<cln::cl_RA>(value)));
	
	else if (!this->is_real()) {  // complex case, handle Q(i):
		const cln::cl_RA r = cln::the<cln::cl_RA>(cln::realpart(value));
		const cln::cl_RA i = cln::the<cln::cl_RA>(cln::imagpart(value));
		if (cln::instanceof(r, cln::cl_I_ring) && cln::instanceof(i, cln::cl_I_ring))
			return numeric(*this);
		if (cln::instanceof(r, cln::cl_I_ring) && cln::instanceof(i, cln::cl_RA_ring))
			return numeric(cln::complex(r*cln::denominator(i), cln::numerator(i)));
		if (cln::instanceof(r, cln::cl_RA_ring) && cln::instanceof(i, cln::cl_I_ring))
			return numeric(cln::complex(cln::numerator(r), i*cln::denominator(r)));
		if (cln::instanceof(r, cln::cl_RA_ring) && cln::instanceof(i, cln::cl_RA_ring)) {
			const cln::cl_I s = cln::lcm(cln::denominator(r), cln::denominator(i));
			return numeric(cln::complex(cln::numerator(r)*(cln::exquo(s,cln::denominator(r))),
			   	   	                    cln::numerator(i)*(cln::exquo(s,cln::denominator(i)))));
		}
	}
	// at least one float encountered
	return numeric(*this);
}


/** Denominator.  Computes the denominator of rational numbers, common integer
 *  denominator of complex if real and imaginary part are both rational numbers
 *  (i.e denom(4/3+5/6*I) == 6), one in all other cases. */
const numeric numeric::denom() const
{
	if (cln::instanceof(value, cln::cl_I_ring))
		return *_num1_p;  // integer case
	
	if (cln::instanceof(value, cln::cl_RA_ring))
		return numeric(cln::denominator(cln::the<cln::cl_RA>(value)));
	
	if (!this->is_real()) {  // complex case, handle Q(i):
		const cln::cl_RA r = cln::the<cln::cl_RA>(cln::realpart(value));
		const cln::cl_RA i = cln::the<cln::cl_RA>(cln::imagpart(value));
		if (cln::instanceof(r, cln::cl_I_ring) && cln::instanceof(i, cln::cl_I_ring))
			return *_num1_p;
		if (cln::instanceof(r, cln::cl_I_ring) && cln::instanceof(i, cln::cl_RA_ring))
			return numeric(cln::denominator(i));
		if (cln::instanceof(r, cln::cl_RA_ring) && cln::instanceof(i, cln::cl_I_ring))
			return numeric(cln::denominator(r));
		if (cln::instanceof(r, cln::cl_RA_ring) && cln::instanceof(i, cln::cl_RA_ring))
			return numeric(cln::lcm(cln::denominator(r), cln::denominator(i)));
	}
	// at least one float encountered
	return *_num1_p;
}


/** Size in binary notation.  For integers, this is the smallest n >= 0 such
 *  that -2^n <= x < 2^n. If x > 0, this is the unique n > 0 such that
 *  2^(n-1) <= x < 2^n.
 *
 *  @return  number of bits (excluding sign) needed to represent that number
 *  in two's complement if it is an integer, 0 otherwise. */    
int numeric::int_length() const
{
	if (cln::instanceof(value, cln::cl_I_ring))
		return cln::integer_length(cln::the<cln::cl_I>(value));
	else
		return 0;
}

//////////
// global constants
//////////

/** Imaginary unit.  This is not a constant but a numeric since we are
 *  natively handing complex numbers anyways, so in each expression containing
 *  an I it is automatically eval'ed away anyhow. */
const numeric I = numeric(cln::complex(cln::cl_I(0),cln::cl_I(1)));


/** Exponential function.
 *
 *  @return  arbitrary precision numerical exp(x). */
const numeric exp(const numeric &x)
{
	return cln::exp(x.to_cl_N());
}


/** Natural logarithm.
 *
 *  @param x complex number
 *  @return  arbitrary precision numerical log(x).
 *  @exception pole_error("log(): logarithmic pole",0) */
const numeric log(const numeric &x)
{
	if (x.is_zero())
		throw pole_error("log(): logarithmic pole",0);
	return cln::log(x.to_cl_N());
}


/** Numeric sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical sin(x). */
const numeric sin(const numeric &x)
{
	return cln::sin(x.to_cl_N());
}


/** Numeric cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical cos(x). */
const numeric cos(const numeric &x)
{
	return cln::cos(x.to_cl_N());
}


/** Numeric tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical tan(x). */
const numeric tan(const numeric &x)
{
	return cln::tan(x.to_cl_N());
}
	

/** Numeric inverse sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical asin(x). */
const numeric asin(const numeric &x)
{
	return cln::asin(x.to_cl_N());
}


/** Numeric inverse cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical acos(x). */
const numeric acos(const numeric &x)
{
	return cln::acos(x.to_cl_N());
}
	

/** Numeric arcustangent.
 *
 *  @param x complex number
 *  @return atan(x)
 *  @exception pole_error("atan(): logarithmic pole",0) if x==I or x==-I. */
const numeric atan(const numeric &x)
{
	if (!x.is_real() &&
	    x.real().is_zero() &&
	    abs(x.imag()).is_equal(*_num1_p))
		throw pole_error("atan(): logarithmic pole",0);
	return cln::atan(x.to_cl_N());
}


/** Numeric arcustangent of two arguments, analytically continued in a suitable way.
 *
 *  @param y complex number
 *  @param x complex number
 *  @return -I*log((x+I*y)/sqrt(x^2+y^2)), which is equal to atan(y/x) if y and
 *    x are both real.
 *  @exception pole_error("atan(): logarithmic pole",0) if y/x==+I or y/x==-I. */
const numeric atan(const numeric &y, const numeric &x)
{
	if (x.is_zero() && y.is_zero())
		return *_num0_p;
	if (x.is_real() && y.is_real())
		return cln::atan(cln::the<cln::cl_R>(x.to_cl_N()),
		                 cln::the<cln::cl_R>(y.to_cl_N()));

	// Compute -I*log((x+I*y)/sqrt(x^2+y^2))
	//      == -I*log((x+I*y)/sqrt((x+I*y)*(x-I*y)))
	// Do not "simplify" this to -I/2*log((x+I*y)/(x-I*y))) or likewise.
	// The branch cuts are easily messed up.
	const cln::cl_N aux_p = x.to_cl_N()+cln::complex(0,1)*y.to_cl_N();
	if (cln::zerop(aux_p)) {
		// x+I*y==0 => y/x==I, so this is a pole (we have x!=0).
		throw pole_error("atan(): logarithmic pole",0);
	}
	const cln::cl_N aux_m = x.to_cl_N()-cln::complex(0,1)*y.to_cl_N();
	if (cln::zerop(aux_m)) {
		// x-I*y==0 => y/x==-I, so this is a pole (we have x!=0).
		throw pole_error("atan(): logarithmic pole",0);
	}
	return cln::complex(0,-1)*cln::log(aux_p/cln::sqrt(aux_p*aux_m));
}


/** Numeric hyperbolic sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical sinh(x). */
const numeric sinh(const numeric &x)
{
	return cln::sinh(x.to_cl_N());
}


/** Numeric hyperbolic cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical cosh(x). */
const numeric cosh(const numeric &x)
{
	return cln::cosh(x.to_cl_N());
}


/** Numeric hyperbolic tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical tanh(x). */
const numeric tanh(const numeric &x)
{
	return cln::tanh(x.to_cl_N());
}
	

/** Numeric inverse hyperbolic sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical asinh(x). */
const numeric asinh(const numeric &x)
{
	return cln::asinh(x.to_cl_N());
}


/** Numeric inverse hyperbolic cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical acosh(x). */
const numeric acosh(const numeric &x)
{
	return cln::acosh(x.to_cl_N());
}


/** Numeric inverse hyperbolic tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical atanh(x). */
const numeric atanh(const numeric &x)
{
	return cln::atanh(x.to_cl_N());
}


/*static cln::cl_N Li2_series(const ::cl_N &x,
                            const ::float_format_t &prec)
{
	// Note: argument must be in the unit circle
	// This is very inefficient unless we have fast floating point Bernoulli
	// numbers implemented!
	cln::cl_N c1 = -cln::log(1-x);
	cln::cl_N c2 = c1;
	// hard-wire the first two Bernoulli numbers
	cln::cl_N acc = c1 - cln::square(c1)/4;
	cln::cl_N aug;
	cln::cl_F pisq = cln::square(cln::cl_pi(prec));  // pi^2
	cln::cl_F piac = cln::cl_float(1, prec);  // accumulator: pi^(2*i)
	unsigned i = 1;
	c1 = cln::square(c1);
	do {
		c2 = c1 * c2;
		piac = piac * pisq;
		aug = c2 * (*(bernoulli(numeric(2*i)).clnptr())) / cln::factorial(2*i+1);
		// aug = c2 * cln::cl_I(i%2 ? 1 : -1) / cln::cl_I(2*i+1) * cln::cl_zeta(2*i, prec) / piac / (cln::cl_I(1)<<(2*i-1));
		acc = acc + aug;
		++i;
	} while (acc != acc+aug);
	return acc;
}*/

/** Numeric evaluation of Dilogarithm within circle of convergence (unit
 *  circle) using a power series. */
static cln::cl_N Li2_series(const cln::cl_N &x,
                            const cln::float_format_t &prec)
{
	// Note: argument must be in the unit circle
	cln::cl_N aug, acc;
	cln::cl_N num = cln::complex(cln::cl_float(1, prec), 0);
	cln::cl_I den = 0;
	unsigned i = 1;
	do {
		num = num * x;
		den = den + i;  // 1, 4, 9, 16, ...
		i += 2;
		aug = num / den;
		acc = acc + aug;
	} while (acc != acc+aug);
	return acc;
}

/** Folds Li2's argument inside a small rectangle to enhance convergence. */
static cln::cl_N Li2_projection(const cln::cl_N &x,
                                const cln::float_format_t &prec)
{
	const cln::cl_R re = cln::realpart(x);
	const cln::cl_R im = cln::imagpart(x);
	if (re > cln::cl_F(".5"))
		// zeta(2) - Li2(1-x) - log(x)*log(1-x)
		return(cln::zeta(2)
		       - Li2_series(1-x, prec)
		       - cln::log(x)*cln::log(1-x));
	if ((re <= 0 && cln::abs(im) > cln::cl_F(".75")) || (re < cln::cl_F("-.5")))
		// -log(1-x)^2 / 2 - Li2(x/(x-1))
		return(- cln::square(cln::log(1-x))/2
		       - Li2_series(x/(x-1), prec));
	if (re > 0 && cln::abs(im) > cln::cl_LF(".75"))
		// Li2(x^2)/2 - Li2(-x)
		return(Li2_projection(cln::square(x), prec)/2
		       - Li2_projection(-x, prec));
	return Li2_series(x, prec);
}

/** Numeric evaluation of Dilogarithm.  The domain is the entire complex plane,
 *  the branch cut lies along the positive real axis, starting at 1 and
 *  continuous with quadrant IV.
 *
 *  @return  arbitrary precision numerical Li2(x). */
const numeric Li2(const numeric &x)
{
	if (x.is_zero())
		return *_num0_p;
	
	// what is the desired float format?
	// first guess: default format
	cln::float_format_t prec = cln::default_float_format;
	const cln::cl_N value = x.to_cl_N();
	// second guess: the argument's format
	if (!x.real().is_rational())
		prec = cln::float_format(cln::the<cln::cl_F>(cln::realpart(value)));
	else if (!x.imag().is_rational())
		prec = cln::float_format(cln::the<cln::cl_F>(cln::imagpart(value)));
	
	if (value==1)  // may cause trouble with log(1-x)
		return cln::zeta(2, prec);
	
	if (cln::abs(value) > 1)
		// -log(-x)^2 / 2 - zeta(2) - Li2(1/x)
		return(- cln::square(cln::log(-value))/2
		       - cln::zeta(2, prec)
		       - Li2_projection(cln::recip(value), prec));
	else
		return Li2_projection(x.to_cl_N(), prec);
}


/** Numeric evaluation of Riemann's Zeta function.  Currently works only for
 *  integer arguments. */
const numeric zeta(const numeric &x)
{
	// A dirty hack to allow for things like zeta(3.0), since CLN currently
	// only knows about integer arguments and zeta(3).evalf() automatically
	// cascades down to zeta(3.0).evalf().  The trick is to rely on 3.0-3
	// being an exact zero for CLN, which can be tested and then we can just
	// pass the number casted to an int:
	if (x.is_real()) {
		const int aux = (int)(cln::double_approx(cln::the<cln::cl_R>(x.to_cl_N())));
		if (cln::zerop(x.to_cl_N()-aux))
			return cln::zeta(aux);
	}
	throw dunno();
}

class lanczos_coeffs
{
	public:
		lanczos_coeffs();
		bool sufficiently_accurate(int digits);
		int get_order() const { return current_vector->size(); }
		cln::cl_N calc_lanczos_A(const cln::cl_N &) const;
	private:
		// coeffs[0] is used in case Digits <= 20.
		// coeffs[1] is used in case Digits <= 50.
		// coeffs[2] is used in case Digits <= 100.
		// coeffs[3] is used in case Digits <= 200.
		static std::vector<cln::cl_N> *coeffs;
		// Pointer to the vector that is currently in use.
		std::vector<cln::cl_N> *current_vector;
};

std::vector<cln::cl_N>* lanczos_coeffs::coeffs = 0;

bool lanczos_coeffs::sufficiently_accurate(int digits)
{	if (digits<=20) {
		current_vector = &(coeffs[0]);
		return true;
	}
	if (digits<=50) {
		current_vector = &(coeffs[1]);
		return true;
	}
	if (digits<=100) {
		current_vector = &(coeffs[2]);
		return true;
	}
	if (digits<=200) {
		current_vector = &(coeffs[3]);
		return true;
	}
	return false;
}

cln::cl_N lanczos_coeffs::calc_lanczos_A(const cln::cl_N &x) const
{
	cln::cl_N A = (*current_vector)[0];
	int size = current_vector->size();
  	for (int i=1; i<size; ++i)
     	A = A + (*current_vector)[i]/(x+cln::cl_I(-1+i));
	return A;
}

// The values in this function have been calculated using the program
// lanczos.cpp in the directory doc/examples. If you want to add more
// digits, be sure to read the comments in that file.
lanczos_coeffs::lanczos_coeffs()
{	if (coeffs)
		return;
	/* Use four different arrays for different accuracies. */
	coeffs = new std::vector<cln::cl_N>[4];
	std::vector<cln::cl_N> coeffs_12(12);
	/* twelve coefficients follow. */
	coeffs_12[0] = "1.000000000000000002194974863102775496587";
	coeffs_12[1] = "133550.502942477423232096703994753698903";
	coeffs_12[2] = "-492930.93529936026920053070245469905582";
	coeffs_12[3] = "741287.473697611642492293025524275986598";
	coeffs_12[4] = "-585097.37760399665198416642641725036094";
	coeffs_12[5] = "260425.270330385275465083772352301818652";
	coeffs_12[6] = "-65413.3533961142651069690504470463782994";
	coeffs_12[7] = "8801.45963508441793636152568413199291892";
	coeffs_12[8] = "-564.805024129362118607692062642312799553";
	coeffs_12[9] = "13.80379833961490898061357227729422691903";
	coeffs_12[10] = "-0.0807817619724537563116612761921260762075";
	coeffs_12[11] = "3.47974801622326717770813986587340515986E-5";
	coeffs[0].swap(coeffs_12);
	std::vector<cln::cl_N> coeffs_30(30);
	/* thirty coefficients follow. */
	coeffs_30[0] = "1.0000000000000000000000000000000000000000000000445658922238202528026977308762";
	coeffs_30[1] = "1.40445649204966682962030786915579421135474600150789821268713805046080310901683E13";
	coeffs_30[2] = "-1.4473384178280338809560100504713144673757322488310852336205875273000116908753E14";
	coeffs_30[3] = "6.9392104219998816400402602197781299548036066538116472480223222192156630720206E14";
	coeffs_30[4] = "-2.05552680548452350127164925238339710431333013110755662640014074226849466382297E15";
	coeffs_30[5] = "4.21346047774975891986783355395961145235696863271597017695734168781011785582523E15";
	coeffs_30[6] = "-6.3439111294220458481092019992445750626799029041090235945435769621790257585491E15";
	coeffs_30[7] = "7.2684029986336427327225410026373012514882246322145965580608264703248155838791E15";
	coeffs_30[8] = "-6.4784969409198000751978874152931803231807770528527455966624850088042561231024E15";
	coeffs_30[9] = "4.5545745239457403086706103662737668418631761744785802123770605916210445083544E15";
	coeffs_30[10] = "-2.54592491966737919409139938046543941491145224466411852277136834553178078105403E15";
	coeffs_30[11] = "1.1356718195163150156198936885250451780214219874255251444701005988134747787666E15";
	coeffs_30[12] = "-4.04275236298036712070700727222520609783336229393218886420197964965371362011123E14";
	coeffs_30[13] = "1.14472757259832757229433124273590647229089622322597383276758880048004748372644E14";
	coeffs_30[14] = "-2.56166271828342920179612184110684658183432315551120625854181503468327037516717E13";
	coeffs_30[15] = "4.4861708254018935131376878973710146069395814469656232761173409397653101421558E12";
	coeffs_30[16] = "-6.0657495816705687896607821799338217335976369800808791959096705890743701166037E11";
	coeffs_30[17] = "6.21975328147406581536747878587069711930541459818297675578654403265380823122363E10";
	coeffs_30[18] = "-4.7255003764027411113501086372508071116675161078057298991208060427341079636661E9";
	coeffs_30[19] = "2.5814613908651936680441351265410235295992556406609945442133129515256889464315E8";
	coeffs_30[20] = "-9752115.5047412418881417732027953903591189993329461844657371497174389592441887";
	coeffs_30[21] = "242056.60372411758318197954509546521913927205056839365620249547101194072057318";
	coeffs_30[22] = "-3686.17673045938850138289555088011327333352145765167200561022138925168680049115";
	coeffs_30[23] = "31.3494924501834034405048975310989414795238339283146314931357877820190435258517";
	coeffs_30[24] = "-0.130254774344853676030752542814176943723937677940441021884132211221409382350105";
	coeffs_30[25] = "2.16625679868432886771581352257834967866602495378408740265571976698475288337338E-4";
	coeffs_30[26] = "-1.05077239977528252603869373455592388508233760416601143477182890107978206726294E-7";
	coeffs_30[27] = "8.5728436055212340846907439451102962820713733082683634385104363203776378266115E-12";
	coeffs_30[28] = "-3.9175430218003196379961975369936752665267219444417121562332986822123821080906E-17";
	coeffs_30[29] = "1.06841715008998384033789050831892757796251622802680860264598247667384268519263E-24";
	coeffs[1].swap(coeffs_30);
	std::vector<cln::cl_N> coeffs_60(60);
	/* sixty coefficients follow. */
	coeffs_60[0] = "1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000007301368866363013444179014835363181183419450549774";
	coeffs_60[1] = "2.13152397525281235754468356918725048606852617746577461250754322057711822075135461598274984226013367948201688447853106595646692682568953E26";
	coeffs_60[2] = "-4.548529924829267669336610112411669181387790087825260737133755173032543313325682598833009521765336124891170163525664509845740222794717604E27";
	coeffs_60[3] = "4.6879437426294973235875133160595324795437824160731608900005486977197800919261614723948577079551305728583507312310069280623018775850412E28";
	coeffs_60[4] = "-3.10861265267020467624457768823845414206135580030123228715133927538323570190367768297139526311161786169387040978744732051184844409191231E29";
	coeffs_60[5] = "1.490599577483981276717037178787147902256911799467742317379590487947009001487476793680630580522955117318124168494382267800788736334308E30";
	coeffs_60[6] = "-5.50755504045738806940255910881807353185463857314393682608295373644157298562106198431098170107741597645409216199785852260920496247655646E30";
	coeffs_60[7] = "1.631668518639067070100242032960081591016027803392225476881353619523143028349554534276268195490790113905102273979193269720381236708853746E31";
	coeffs_60[8] = "-3.9823057865511431381368541930378720290638930941334849821428293955264049587073723565727061718251925950255036781219414607001763225298119E31";
	coeffs_60[9] = "8.16425963140638737297557821827674142140347732117757126331775708561852858085860735359056658172512163756926693444882201094206795155146202E31";
	coeffs_60[10] = "-1.426548236351667330492229413193359354309705120770113917370333660827270957172393778178051742077714657388432785747112574456061555034588373E32";
	coeffs_60[11] = "2.14821861694536170414714365485614715949416083667308573285807894910742621740039595483105992136915471547998283891842897000924199509164799E32";
	coeffs_60[12] = "-2.81233281290021706519566203146379395136352592819625378308636458418501787286411189089807465993150834399778687427813779950602826375635436E32";
	coeffs_60[13] = "3.222783358826786224404373038021509245352188734386849874296356404770508945395436142634892645963851510893216093037595555902121365717716154E32";
	coeffs_60[14] = "-3.250409075716999887328836263791911196138647661969351655925350981785153422033954649154242209471752219326556302767677017396179477496948985E32";
	coeffs_60[15] = "2.897783210826628399578158893643627107049805015801395657097255344786041806868455726759715576609013221857885740543509045196763816109465777E32";
	coeffs_60[16] = "-2.29136919195969647663887561122314618826917230275433296293059354280077561407373070937197721317435316121212106870152659174216557412788874E32";
	coeffs_60[17] = "1.611288006928200619663496306945576194382628760891807800193737346171844871295031418730500946186238469256168610033434708290528870722514911E32";
	coeffs_60[18] = "-1.009632466053186015034182792930705530447465885425278324598880797572411588461783484686932989855033967294215840157892487264656571258327313E32";
	coeffs_60[19] = "5.64520651042784179741815642438421132518008517154942873706221206276337451930555926854271086501686252334516011905237101877044320182980053E31";
	coeffs_60[20] = "-2.81912877441595327683492797147781153304080114512116755424671954256427789550109614317215500473322621746416096887803928883800132453510579E31";
	coeffs_60[21] = "1.257934257434294354026338893625531254891110662111965279263894740714811495074726866375858553579650295684850594211744093582249745250079168E31";
	coeffs_60[22] = "-5.01544407232599962845688086323662774702854661522104499328570796808858930542190600193190967249971520736397504227594619670310759235566195E30";
	coeffs_60[23] = "1.786035425040937365122699272239542501767986628253845452136132211710520249195280548478081559036323184490150479070929923213045153333111476E30";
	coeffs_60[24] = "-5.67605430104368150038863866362066081946938075036837029856903803768657069745962581310398542442108872722631658677177822712376500859930109E29";
	coeffs_60[25] = "1.607878222558573982505999018371559631909289246981490321219650132406126936263403946310818841465409950661433241956831540547593847161412447E29";
	coeffs_60[26] = "-4.05332042374309456146169816144083508836132423024788116321074411679252452773181941601763924562378611113519038766273534176937279867894066E28";
	coeffs_60[27] = "9.07493596543985672039002802030098143847503854224661484396413496012780904911929710460264147600378604646912175235271954302119768907744722E27";
	coeffs_60[28] = "-1.800074018924350353143489874038038169034914082090587278672411654146678304871125651069902339241049552886098125667720181441150399048551683E27";
	coeffs_60[29] = "3.154250688078046681602499411296013099183808016176992164829953752437167774310360166977972581670851790753785195101324694758021403186162394E26";
	coeffs_60[30] = "-4.86629244083379932983782216256143990390210226006560452979433243294026128577640975980482675864760717747936401374948595060083674140963469E25";
	coeffs_60[31] = "6.58428611248406176613133080039790689602908099995907522692286902207707012485115422092589779128693214784991500936878932461139361901566087E24";
	coeffs_60[32] = "-7.77846893445970039116628280774361378296946997639645747353868461156972352366479641995295874152354776734003001337605345817120316052066992E23";
	coeffs_60[33] = "7.98268735994772082084918485121285571015813651374688487489679943603727447378945977989630573952891101472578977333720105112837324185659362E22";
	coeffs_60[34] = "-7.07562692971089746095546542541499489835693326760069291570193808615779224025348460132750549389189539682228913778397783434269420284483726E21";
	coeffs_60[35] = "5.381346729881846847476909845563262674288431852755093265786345982700437823098162630059919716651136095720390719236493773958116646152386075E20";
	coeffs_60[36] = "-3.4856856542678356876484367392130359114150104987588151214926676834365219571876912071608359944324610844909103855562977795837329347647911E19";
	coeffs_60[37] = "1.90665542883474657677037950113781854248329048412482665873254624417996252139138481002200079466749149325431679310476862249520001277129217E18";
	coeffs_60[38] = "-8.72254994006151131395107200045641306281165826830744222866994799005490857259177347821280095689079457417603257537321939951004603693393316E16";
	coeffs_60[39] = "3.30066663941625244322555483012774856710545517350986120571194216206848716066355962922968824538055042855044917677713272771363157100391997E15";
	coeffs_60[40] = "-1.020092089391030771746960980075254826475625668908623135552682999358854102567810002206013823466362488147261886160954607897574298699485318E14";
	coeffs_60[41] = "2.537518136375035057088980117582986067754938584307761188810498418760131416720976321039509027979006220650166651208980823946300429957067604E12";
	coeffs_60[42] = "-4.99523339577986301543863423322168947825482352498610406809585164155176248614834684219539096936869521198401912030883142734471627752449382E10";
	coeffs_60[43] = "7.62961024898383965152735310352890448678585029645218309944823403624458716639413808284778269959424212699922000610764015063766429510499158E8";
	coeffs_60[44] = "-8834336.1370238009649936481782352367054397712953420330251745022286767420934395739052638862442455545176778475848478708230456099596423988";
	coeffs_60[45] = "75445.9196169409678879362111492280315111800786619928588067631801224813888137547544321383450353324917130013984795690223150786036557545929";
	coeffs_60[46] = "-459.8458738886001056822131294892698769439281099450630714273592488999986769567563218319365007529495798105783705491469742412340762305916056";
	coeffs_60[47] = "1.922366163948404706136462977961544621491268971185908661903800938507393909575693892375103171073678191394626251633433930639174604982075991";
	coeffs_60[48] = "-0.00524987734300376305383172698735851896799115189212445098242699916121836353753886238290792298378658233479210271064792489583846726184351881";
	coeffs_60[49] = "8.81521840386771771843311455937479573971716020932982441671173279504850522350287085310420429874536637110755391716691475171030099411021337E-6";
	coeffs_60[50] = "-8.42883518072336499031504944519862331274440110738275125460829656821173301216150526266773841539372995424665091651911614576906895281293397E-9";
	coeffs_60[51] = "4.1559932977982056953309753711587342647729282359841592558743510304569204546713517319749817560490538963802716194154620384631597656968764E-12";
	coeffs_60[52] = "-9.26494376646923216540342478135986593801117330292329759013854851055518195892306285985326338987592590319793280515888731024676428929933443E-16";
	coeffs_60[53] = "7.80165274836868312019654872701978288745672229459298320116385383568401529728308916875595120085091565550085090877341856355815270191309086E-20";
	coeffs_60[54] = "-1.922049272463411538721456378153955404697617250978865956250065913541261535132290272529565880980548519758359440057376306817458561627984943E-24";
	coeffs_60[55] = "9.46189821976955264154519811789356895736753858729897267240554901027053652869864043679401817030067356960879571432881603836052222728024736E-30";
	coeffs_60[56] = "-5.06814507370603015985813829025522226614719112357562650414521252967497371724973383019436312018485582224796590023220166954083973156538672E-36";
	coeffs_60[57] = "1.022249951013180267209479446016461291488484443236553319305574600271584296178678167457933405768832443689762998392188667506451117069946568E-43";
	coeffs_60[58] = "-1.158776990252157075591666544736990249102708476419363164106801472497162421792350234416969073422311477683246469337273059290064112071625785E-47";
	coeffs_60[59] = "4.27222387142756413870104074160770434521893587460314314301300261552300727494374933435001642531897059406263033431558827297492879960920275E-49";
	coeffs[2].swap(coeffs_60);
	std::vector<cln::cl_N> coeffs_120(120);
	/* 120 coefficients follow. */
	coeffs_120[0] = "1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000060166025676976656004344991957470171590616719251813003320122316373430327091055571";
	coeffs_120[1] = "3.4497260317073952007403696383770947678893302981614719279265682622766639811173298171511730607823612517530376844024218507032522459279662180470113961839690189982241536061314319614353993672315096520499373373015802582693649149063603309572777186560148513524E52";
	coeffs_120[2] = "-1.4975581565000729527538170857594663742319328925831469933998274880997450758924704742659571258591716460336677591345828722528085692201176737000527729600671680178988361119859420301844184208079614468449296788394212801103162564922199859549237082372776667464E54";
	coeffs_120[3] = "3.1957762163065481328529158845807843312720427291703934903666695190945338610786360201875291048323381336567812569171891600400186742244091402566230953251621720778096033490814848238417212345597975915378369497445590090951446115848410773972658485451963575288E55";
	coeffs_120[4] = "-4.4689623509319752841609439083871154399631153121231062689347162975834499076693093642474289117173045421812089871506249999929076992135798925381959196225961791389783472385803138226317976820364502651110639008585046458007356178875618627927171581950486124233E56";
	coeffs_120[5] = "4.606068718424276543329442566011849623375399823565351941825685310847310447457609082356012685588953435307896055516214072529445026693975872604267789672469025113562486157850515006504573881812473997762948360804814769118883992998548055557441646946685125118E57";
	coeffs_120[6] = "-3.7314461146854666499272326592212099391213696621869706562566612605818861385928266960370453310708465394226398321257508947092784006446784523328681347046673172481746936234783770854350210504707173921547794426833735429199925024679815789545465854297845328325E58";
	coeffs_120[7] = "2.474425401670711256989398808079221298913654027234786607507813220440957186918973475366048940039541074278444160674001228864321389049663140487504402096319272526201782217412803784224929141788255724630940381342478088455751340159338461174261577243566175687E59";
	coeffs_120[8] = "-1.3811875718622847750042362590249762290599823842851465148429257970104907280458901604054390293828410620002370526629527048636126473391278330353375163563724888073254512227198849135923692811222561965740181944727170495185714496890490479692693474125883791901E60";
	coeffs_120[9] = "6.623089858532754482582703479109160446021743439335073883710993620625687271109284320721410901325182604938578905712329203551531862389936804947105415805829404869727743706364603519433193421234231031076682156125442577335383798263985569601899041876776866622E60";
	coeffs_120[10] = "-2.7709515004299938864490083840820063124223529009388282231525445615826433364331567602934962481829061542793349831611106716261513624279121506887680318284535361848032886450351898892264386237450622827397559067350672965967202437971930333676917000390477963866E61";
	coeffs_120[11] = "1.02386112293172223921263435003659366453292875147351461165091656394534393086780717052422266565203902889367201592668259202439166666819852985689989767402099479793087277263747942659943270101657408462079787397068550734516045511611701546009078868077038808757E62";
	coeffs_120[12] = "-3.3740197731917655541744976218513993073761175468772389726802124778433432226803314067431898210976006853342921093194297198044021414900546886804610561082663076825192459864843102368108908666053756409152492134638014803233805912009476407113691438596300794146E62";
	coeffs_120[13] = "9.996217786487670355655374796561399578704294298563457268841140703036898520360123177193155340144551120260016445533739357030180277693840431766824113840895797510199955331557143980793267795747200088810293047731873192410526786931879590684673414288653913515E62";
	coeffs_120[14] = "-2.6804750990199908441443350402311488850281543531194918304012545530803283220192092419107511475988099394746800512008906823331244710178292896561401750166818497729239682879419868799442954945496319510685448344062610897698253544876306888341881254056091234759E63";
	coeffs_120[15] = "6.5422964482833531603879610057815197301035372862466791995455246163778529556707384145891730234453157337303612060344197138180893720879196243783337539071284141345021864817147590781643393947019750353147151780290464319306645652085743359080495090595200531486E63";
	coeffs_120[16] = "-1.4604487304348366496825146715570516556564950771546738885215899741781982964860978993963272314092830563320794184042847908967120542212316261920409301852237223308467032419706968676861616456179895880956772385510853673982424825597152850339588189159102980666E64";
	coeffs_120[17] = "2.9942297466313630831467808691292548682230644559492580161942357031681185068971393754352871129412787966878287389513657398203481589163498279625760316093736277896138061249076616695157422053188087353540756151586375196486093987258640269607104978950906670704E64";
	coeffs_120[18] = "-5.658399283776588293772725313973093187743120982052603944865098913586526167668102733207163739469977271584007101869254711458133873627143366757941713180350370056955237604551850024423889291598422917971467957836705917204903687959901869098540153925178692732E64";
	coeffs_120[19] = "9.887368951584101622633538892976576123080629424367489037686110916264512731398396560326756128205833849608930564615629875435100785011872254223744155330328477703592008501954532369042429700051416733748454350165515933757314793533786385104271308839525639768E64";
	coeffs_120[20] = "-1.6019504550228766725078575508839635707919311420327486864048705642201106895239857903763208049376932672160478820626774934879424498715258948985194011690204294886396827446040036506699933786721588971678753877371518212675519147446728054067639530675249526082E65";
	coeffs_120[21] = "2.4124568469636899540706437405441629738413207418758399778576327598069435452295650039157974716832514441625728576753250737726840109004878753294786785674578138926529507088264657400701828947949531197915861820274684954206665488761473274445827472596875582911E65";
	coeffs_120[22] = "-3.3841653726400000079488483558717068873181168418395106876260246491163166726612427450773591871178866824643300679819366574162583413250423974373322308130319007820863363304629451933781204964221002853140392226489420463827400812929748772154909106349410663293E65";
	coeffs_120[23] = "4.4305670380812288478773282114598811227924298131011853412998479811262358077680067168455361591598296346480072528806092976336961470360354620203822421524751468329936930212919915114854135818230382164555078957880154875221176513434392525189922941290050575762E65";
	coeffs_120[24] = "-5.4228176962574428947233160003094662570284359565811627941401342797491445636152854865132166939274138115146035207618348708829039395974942115203986578386666664945394109693178927438991059414217518334491360514633536224841961444935232548483014691997071543828E65";
	coeffs_120[25] = "6.214554789078092267222051275213928685756510105900211846145778269883351640710249139978059486185007208670776100912863866582278800642692097830681092656540813877576256048148229340562594504915197956922387464825593941922429396202734006609196778697870436014E65";
	coeffs_120[26] = "-6.6773485088895517986512141063848395783979405189075416643094283756118912554557672721632998501682483143868731647507940026369035991063923616298815637819145806214374157182512600196214559297579802178103007615921637577873304407436850546650711237281572008424E65";
	coeffs_120[27] = "6.734925330317694704469314845373778479111077864464012553672292377883525864326847400954413754291163739900219432201437895152976917857427306427115230048061308424221525123820493252697918698598513232640014129066982507718245232516657455821629338155744427538E65";
	coeffs_120[28] = "-6.383582878496429871173501676061533991181960023885889537277705274319508246322757005217436814481703326467002683699047193244918123789600842413060331898515872574523803039779899326755393070345055586059441271293717500426377884349137309244757708993087958455E65";
	coeffs_120[29] = "5.6913529405959275511022780614007027176288843526260372650173869440228336395668389555081751187360483397341349300975285817498083216487282169140596290796279875175764991375447348355187090404486257481827615256024271536396461908482904537799521891879785332367E65";
	coeffs_120[30] = "-4.776996734587211249165031248400648409423153869394988746115380756083311986805300361459383722536698540926452976310737678416019979202990255666249869917768885659350216400547190883549730059461513588008706974008085270389354525600694962352952715682056375518E65";
	coeffs_120[31] = "3.7775367524287124255443145064623569295746034916892464094281613465046063954544055573214473155479196552309207209647540614474216985097792266203411723082566949978062697757983354600199199984244302856099811940389910544756210676400851240882142140969864764468E65";
	coeffs_120[32] = "-2.8161966171919236962021901287860232075259781334554793534017516884995332700369401674058517414969240048359891934178343992080557338603528540157030635217682829894098359736903943078409166055640608627322968554856315650475005062493399450913753277478547118352E65";
	coeffs_120[33] = "1.9804651678733327456903212258413521470612733719543558365536494344764973229749132899499862883369665827727506916597326744330471802598610837032598656197205238983585794266213317465548361566327435762497208877015986690267754534342053368396181078097467171858E65";
	coeffs_120[34] = "-1.3144275813663231527166312401997093907605894997476799416306355417933431514642211250592825223377757973148122542735038736133300194844844655961425683877005418926364412294006123974642296395931311307760050290069031276972832755406161248577410950671224318855E65";
	coeffs_120[35] = "8.2367260670024829522614096155108151082106397954565823313893008773930966293786646885943761866773022391428854862805955553810619924412431932999726399857050871862122529700098570542876369425991842818202826823540112018849926644955200888291063471724203391548E64";
	coeffs_120[36] = "-4.8749964750377069822933994525197085013480654713783888755556109773660249389776804499013517227967180500633060271953473316017147397601291325922904139209860429881054757911243087427393920494271315804033914011087815785282473032714919188637172020633929566123E64";
	coeffs_120[37] = "2.7259664773094932979328467102942769029907299417406744864696200699394594868759231280169149208728483197299648608091313447896342349454038879581019820193316159535211365363553004387852005780736869678460092714636910972426808304270369152189989142121207224142E64";
	coeffs_120[38] = "-1.440426226855027726783521340050349148103881707415523724377763633849488875095817796257895327883428230885349760692732068174527147156893314818583058727424251827006457849321094911262818557954829248070170426870959233263267490276774734065709978749400927185E64";
	coeffs_120[39] = "7.193790858249547212173205531149034887209275529426061411129294234841122474820371873361100215884757249851960370114629943083807936135915003201800204713978377250292881453568756354858194614039311248345228434431020394729104593125888325843724239404594830488E63";
	coeffs_120[40] = "-3.3960029336234301755324970935705944032408435186630159101426062821929524761770439420961993430248258136340087498829339209014794230274407979103789924433683527009234592433480445831820377517333956042612961562022604325181492952329031432513768020816986814393E63";
	coeffs_120[41] = "1.5154618904112106565112797443687014834429200069480460967081898435635890576815349145926430052596468033907024005478559584915319911380449387176530845634833237204659108290330613043367085829373476690728522550189678729181372902816898536141595215616716630939E63";
	coeffs_120[42] = "-6.3927647843464050458917092484911245813170740434503951669888756878206365814594631676413018245438405308353724023007754523096143775098898268650326908751515418318201372985246418468844138298345777180517875695389655616000832495210812684049030674085212697428E62";
	coeffs_120[43] = "2.5490394366379355452002449693074954071810215414182359403355645652443600688717811337587901850157210686351097591461582890354662732336749618027675479531031836144519267481752770036252137747675754903974915999567019837855523058289177148692481402871253211324E62";
	coeffs_120[44] = "-9.606466879185328464666445215840505657671157752044466089989040292763536710311599947887918708456526669882072519263973105599580140713596301388561639705589314111762600854460059589939760935803484446888352368360433606245369171819922425771642408570388554052E61";
	coeffs_120[45] = "3.4212392804723358445152430359637323789304939688937873921904941412927756295848104328630952153624979489607759834359194243032109828811134607612715016533909375981353098879969472700079242226099049323998286020341979178782935852542355220403299144612362244738E61";
	coeffs_120[46] = "-1.15119134701605919461057899755821946453925102458815313053351247263978303346790555715641513756644038607667203392289423834966320935498856390723555744530789850290369071103208529608463210398231590077340268751005311531529356083188150256469829678612245577446E61";
	coeffs_120[47] = "3.6588622583411432033523084711047679684233731128914152509273818448610176621874654252431411048902598388935105893946323641003087000410095802098177375492833543391040706755511234323104846636415419597151008153829618275459606044459923718022154035121167198784E60";
	coeffs_120[48] = "-1.0981148514467449476248066871827754422009180048705085132882492434176164929454140182025449006310206725429473330511884213470600326782740663313311256352613244044500057688932314549669435095761340307817735687643806167483576999980691227831561891265615422486E60";
	coeffs_120[49] = "3.110998209448997739767747906196101611409160829345058138064861244336130082424927251851805875584197897229644157110035272012393338413235208343342708685139492629786435072305986349067452739209758702078026647999828517440754895711519542954337931090643534216E59";
	coeffs_120[50] = "-8.3162485922574890007748232799240657004521608654422032389269811102140449056333167761296051794842882201869698963586030628312914066893199727852512779320175952962772072653493447297721128265231294406156925752496310087025926300388984242024436858845487466277E58";
	coeffs_120[51] = "2.0966904516945699848169820408710416999765756367767199815424586610234585829069218729220161654233351574517459523275756901094737085187558904179251813051891939079067686519817858153690134828671544815635956527611986498479411756457222935682849773436423295467E58";
	coeffs_120[52] = "-4.983121305881207125553776640558094509942884568949257704810973508397697839859902664482541160531856121365759763455699578413261749913567077796586919935391984240753355552646184306812426079133011894826183873855851966310877619118554510972675999316631346679E57";
	coeffs_120[53] = "1.1158023601951707374356047495258406415892974604387009613173591921419195864040428221070481312383179580486787822935456571355463718115785982888531393271665510645725283439572279946304699780331972095822869500426555507626639723865965516308476400920600382357E57";
	coeffs_120[54] = "-2.3524850615012075127499506758220926725372558166170912192116695445007095502575329450463479860779122789467638956004572617263549199692255055063165454868102165975951768676031140009643202074220557325155838768661030361538572755082660730808847591840060467064E56";
	coeffs_120[55] = "4.6669318431895615057208826641721251136909284138581355667925884903657855204100373961676117747969449100495897986226609480142908763981931305129946569690612924941456739524153327260627627771254850382983581593260532259539447965597396206625726656509884058042E55";
	coeffs_120[56] = "-8.7053773217442419007560462613131691749734845382618514999712446313788486289774350240165530159591402631439776213579542026449818009956904779042347595401565525081115611496250192338958392965746523979241969677734430475813057146574920495171984815351708574336E54";
	coeffs_120[57] = "1.5256552489620511464542280446639568546874380361953025589702692266626310669215652044048704882910412155084167930513006634430352568411276836880182348033924636960897794333644980768878022821035659978039286230061734024129667272393315169114199838321062607299E54";
	coeffs_120[58] = "-2.5099934505534008439782195609383796207770494575364994376922414269548303512602084430128307108303305643530918354709126474742035537827601791192999467996479881350277448357927640707861695639576629921988481117017137420422963638277868648516492581097660522547E53";
	coeffs_120[59] = "3.872963359882179682964169603201046384616694634651871844057456079738892419308420856725974686574980381399016464501318163662938118593626674643538005780375691959391996340141057698193381380484420715733863044826589570328349973407598034428591146829028071358E52";
	coeffs_120[60] = "-5.59947633823301408044455223877913062308847941596689956112764416031828413291312481723036534632655608672535030921469531903033364444816678754679807809159478411100820014592865068932440734964265842594875758737421026093110624848762070026616564150314951394E51";
	coeffs_120[61] = "7.57762861280525531438216991274899157834431478755285945898172885086150762425529113816148806028462888396660067975773261101497666568988246606837690320098870044112671149076084444095163491848634465373822951831018725769263871497616640732007420499659069842E50";
	coeffs_120[62] = "-9.587786106526273406187878833167940811862067040706459726637556599860244751467528905534431960251166924163661188573831350928972391892492380823531476387272791432306808700507685765850397294118719242350333451452137838374120658600691461454898577711260078952E49";
	coeffs_120[63] = "1.13288726401696728230264357306938076698155303500407071418573081766541065136778223998897791839613776442037036668986628122296219518360439574147622758002647495909592177914657175019781723803408732148262293125845657503039410078589916085532057725749397276232E49";
	coeffs_120[64] = "-1.24849787223197441956280303618704887038709792250544105638342097080498907831514597860418910331910245753340059089147824955071899315894649696314820492532126554883819507650973976145456660786429117569053901704116877128391672511345177517877672824534972448216E48";
	coeffs_120[65] = "1.2815463720972693091316233381473056495608681859925407504190742949467232967966271661733907550222983737930524555721493736920130260377888287772008209963158064973076933575966719577456540496444474944074979736374259087350416613616719928507635667369740203319E47";
	coeffs_120[66] = "-1.2234887340201843394744986892310393596065877342193196880417674427168862926389642850813687099959036354499094230765541977493433449153438766822382486040215211159359175689369230076522107734270943423777076523650345103234411047700646432924770659676420158487E46";
	coeffs_120[67] = "1.0847187881607033339631651118075716564835185723270640503055198532318419482330026641941088359447807553514405522074008969583213861070993661224871455023365601323302778638456843760403418046238489404394483720438784739822580385277055304353975028280477740796E45";
	coeffs_120[68] = "-8.9160881476675795743767277986448579964735858351472748620623279571408606135698760493224031735408212513500922230670883171668702983221921543376953865813604783695111225412173880768170509738290662806468458720236121755965944855709552219268353813402612336565E43";
	coeffs_120[69] = "6.782864920104031936272293608616215844503387641476821968620772153274069873138756405621471099960069602613619793775294358177761533027360002770186566164041138064221354961783144649476276625776241973967317262115970868665380343599565811072109785000646703404E42";
	coeffs_120[70] = "-4.7667808452660756441368384708874451089976319738852731080495062883240643961463680300964077232336439626019128672679703771884184482488932861160134911816225569323838390204451496983578077563176966732010513231048738892639707790407292070646798259086924770995E41";
	coeffs_120[71] = "3.0885057140860079424719232591765602418793465632939298397987628606701994268384966881159469651774584648643122830739130127593326652998108850492039117928976417052691273804304806596509726701594300563830431015215234640024338277573401498998072908815285293868E40";
	coeffs_120[72] = "-1.8410405906573614531857309495652487774337134256805076777639383854080936219680656594060736479739035202182601529001321266214227848431889644620036213870966329509961114940541333851155401637197303308322414678191211465563854205816313387785764908216851396633E39";
	coeffs_120[73] = "1.0073694433024942271325653907485159683302928496826793112696958500366488338508211620934892875328717073528902110227362794694820010124321343709182901273795782541866547318841893692957109947576483162095037812781379193423759617638948859880051822460818418552E38";
	coeffs_120[74] = "-5.0475051506252944853315611134428802424958512917967945464108691542854207821486654807141339210375899950551724141366521361887864357385178212628348794663127149312605456165451981719848656127310229221238908657530297751682848475855876378576874607521597136906E36";
	coeffs_120[75] = "2.3099766115359817610656986443137072041797751710805647712896098246833051023271876304983288225638204962631413469467959017768113430777226924099787875749611560913177631681394153889301715579572842026181746028117354815826836594637709952294015960031772162547E35";
	coeffs_120[76] = "-9.629053850440590569772960665435833408449876392175761493622541259322053209458881628458334353756739601360772251654643632187697620334088992038575944303101187678397564511853344433267011583960451100374611538881978045643233876974513962362084978095067025623E33";
	coeffs_120[77] = "3.6452126546120530579393646694066971671091434168707822859890104373691687449831950255953317231572802167174179528347370588567969602221261721708890001616085516755796796282628169745443137768549800602834096924025507345446292715781107949529692160434800323E32";
	coeffs_120[78] = "-1.2492564030201607643388368733220662634846470405464496879151879822123866671204541555507638492613046717628358162773937737774832271305618491107140304474323049182605167775847584622690299098207979849043605983558768056117581593008210986863088433891075743152E31";
	coeffs_120[79] = "3.8627447638297686357472526935538070834588578920414538227245516723308987020816841052950727259618753144711425856434270832495754300189881199851254605718213699755258867641301730599979474865704144160112269948588154919128986989885090481959424806312935273075E29";
	coeffs_120[80] = "-1.0736758703963497284148841547397192249226725101007524773889805877171959717011395181953504058502607435217886087332761920207901621377557079619638699346496468750455986591040017334237734940082333954589067611955107878899677189289648293223359861027746438121E28";
	coeffs_120[81] = "2.6722714785740082059347577649909834926335247252399259683264830680945466475595847553753509546415283809619181144796536494882020159787371993099998263815645014317923922311421330376008111312767167437401741178863083976628261471599264811824656877164988491393E26";
	coeffs_120[82] = "-5.9304047185329750657632568788530498935629656326502947505210292278638825286675833282579834326765999907183142489791905921257123760969245535649745876992946512033156167841406724363867902645010435996961270021857807247440211477908060243655541266857227638988E24";
	coeffs_120[83] = "1.16817022089143274700208191285335154155497013626172270535715899131321522799010543339535307798264602677955894930046454353008462671803498794203612585729705145312299224155123919877760274781582850868001155383467754608529345730226972329454404720862870618607E23";
	coeffs_120[84] = "-2.03239515657536501213472165328009690017090356606547792466197690386716728380893226886179282271040418637806139515373566132123131620086873213475424131345589653019635327048678766191769576650893957440830876852296666120473866301097954633389040518870395767125E21";
	coeffs_120[85] = "3.1065334503269182605978912331263087603258864771943471481540265718169544724355602987297631515907391374943512439350265433478241465606056187134785807375293801936399644663199667496663518171930757047012102683120173353568660795955174938680248863153863947508E19";
	coeffs_120[86] = "-4.1476244154347831048636005592317388215032295704489937704602030038303705695463546496640638505584602502764898113504560236629804442607426019604639559875021291459916615723777004493344143132459204229291886967479716413925814352313734234340863490128872380307E17";
	coeffs_120[87] = "4.8067293487250079673131214670887682215073707729621636364424152483295071605326220176372385638491275365750175037404843071051780212494354459897540110089573898336327006157766256896984455454193271731091632286742192439925748114360605084629432813597189767538E15";
	coeffs_120[88] = "-4.8023544548381246628003457039588616467438691159189277447469028024236284353593054364114519649309416187375157096932150251663679454372678125518452171003992957433311257042292636706448339781439297178835786059318810522834929923770539615271536113963729385909E13";
	coeffs_120[89] = "4.1055087514683476865343055835875083237542317413651906253058979029083965525058905726360233143503628224856307545474786181299719957472120906835233967660557875100202077212004953379299507351564181758434304881046845705855303854083493519588411179065109026834E11";
	coeffs_120[90] = "-2.9787503393847675871205038539267895335240592213878943742323972872214441728681744433089698206110260166068266926018988659692353298939109421567999207730700359726920482465669373553804927535369930188390246988033893916611435406224816632683980860607732310186E9";
	coeffs_120[91] = "1.8178328110729629877907010659834277046059726898311908447099830056830012488194646687474150289147446390570639168063598563291822008033517936194534129929881015025633519502485415000390171249019651579295905194415531994026553693578406432674734610095421683863E7";
	coeffs_120[92] = "-92391.136314434380495997449781381513978328604842061708454699991154771188446049720445502194923435235472458378926242100033122111143321209059959788378033220861638093951546784186137626553022963832352496255851690092415165826965388502958309163995296640164754";
	coeffs_120[93] = "386.82763074890451546182061419449593717951707520472938425276820204065379182568600735469831672149785863654956632602671563997131280046154927653332261114114005498875447205079045401364007035880825957300757663780618819785476980699579657587509130753204519233";
	coeffs_120[94] = "-1.3181204292571874302358432444324779303744749959754136019600954846045028319805074783759764870805734807693739252625657350494147444011046941331047057337345953605042408524072436811726898109072388160378243068564382623631658424851676817690976343859083960324";
	coeffs_120[95] = "0.003606538673252695455085947121496196507159591230095595764694813152630524319596509155920374890595867709349176662036024214476302717902368680224618116411588086562230407996267622244422187853090635901906175373997993725355114393033631058067900506212434600015";
	coeffs_120[96] = "-7.805244503909439374422205381130511738566245024242591464192744568789876715121004646510755612128565674260161510215430132815223049297785205382643947556846567064565241387424696940674258789227398935846571768027456535982674711768030751512030174841314425949E-6";
	coeffs_120[97] = "1.31373705470989377112938364152965446631228819123896570245455699237549295870321627234421140232628798373711221392827979836922621437205363811871692678679625916100572037589291239046725228767017131155814257944742981208252138821140381478767814046301821211856E-8";
	coeffs_120[98] = "-1.6872873094408224472617181717534409090015431593544429529131126514352910895332010213914243717484771690790552077128803350550170014347729272790464826195676369023970955260051387240496705602732313607409271794413329062030561818907163134089683283286623809325E-11";
	coeffs_120[99] = "1.6183083251905685095057354853863188515437903228178486856957070037813756492593759658405336450433607296873747595037080703825755020175480385843762609522889527239577435110258291566585028919336090916225831079571865425410181260759913688103716786795647286451E-14";
	coeffs_120[100] = "-1.13097359411474028225398794102354853670936316496817819635688647804142428962171772690717075128208102537660772310780986623575005236651312181907812813813504742701120603881086064664411899253566047514905519888629604717647221817372977488600336785871295304013E-17";
	coeffs_120[101] = "5.599216369109121957949255319730053610385733330502739423509794477602247233276045188197007198417289907263120960704056657544648432653622931077692740961599655386871075693202473992087883485704436336279135221721374640982826144708808646466699352755417123702E-21";
	coeffs_120[102] = "-1.9009180102993021108185348502624676395148544369474718879637745630712451378711342634099259114111847962752555305470572286326367888004493816251811794947276966269738750207359305252041104539066278002044545942171476984766923991983055271262414217352967659228E-24";
	coeffs_120[103] = "4.261262509940940316499754264112111685174274727656165126333137554124192224955656564229887938745508952447664695831728428607673797269945824475565104978593072684829487175697371245288754204324544164474840153141042852153497051337607734150135978754952561336E-28";
	coeffs_120[104] = "-6.033854291373449912236926137860325602686312455380825767485673949251953414778800668020214699151728472172651816317924130614791108454134597377848088327850505473503152696524861086193124979489104732214189466703901268332265826882296309653009237279831825243E-32";
	coeffs_120[105] = "5.1208402745272379096703574714836785944518835939702823617280147111145234914591060871138496110227453241036619229980622243972303295470574470937679143516006222494480144845809123492603651773613707216680534850900104861326332900592715684757980394834998321888E-36";
	coeffs_120[106] = "-2.4463535717946588550832618025289907099586319384566637643650142186828541109926588999585266911960640972919441499109750654299062004147686492034166034659422424525984094382368955916181276646903453872999065929058429821759475215620044891133652431220664754175E-40";
	coeffs_120[107] = "6.0973480699773886324239008989591793773608942051497498591908583910660358857815864266160341286217871697703362816166340947142517661604423899536979689047275448159991318658879804351288744125363072102852651926942302209139318098544348348564409845011546432615E-45";
	coeffs_120[108] = "-7.2234185761285078775026471720270426097727212523472472797635230392183067756271499246654638332288950167477129840028892565652782123508855602380279653475510712205780583313834027906297063690370430285856541927759405826980856379432703473274890527421175151858E-50";
	coeffs_120[109] = "3.6217112680215791206171182969894344487335819731880124290544082848140757826983885738735436324684863867140575000400288923606439193119990961489053513339202655922248092157737577138929144240507796562250602457839068582279379672722261563501188150876583184441E-55";
	coeffs_120[110] = "-6.6329300032795486066608594142675837603786558782159646987663521197523704085781830169369726460621246948945196657495305819768951424025780824076252490918306538895670861455244641773606980519824591785816943621538721352987553804824051051144609050417497894495E-61";
	coeffs_120[111] = "3.6664720904335295532012711597888717227860988776477301054518326674835421172405060906940404374163713097964932859351917152390238690399278248344863365606468942320103392909602843987855082225592776850615943708151738327210634139824601616072015258461809772448E-67";
	coeffs_120[112] = "-4.7466013179695826928232672846686064011594588664906398407027593213652099998530859940288723349213099851532139911079905393494419637612780994270110734378146177806681489226896952731800026849872070824592339117757940119304241732812925979963178130104280115315E-74";
	coeffs_120[113] = "1.0163707785221910939390789816391472677729665860532352695801597334766068288835382195560328979864550624486740471947632369344045378626680607890520366137741785540226552923584183986350590955499329375427326072319268396685478606934920507703868118038891818762E-81";
	coeffs_120[114] = "-3.4814151260242800905467399051937942442621710748397374123807284826536707678408888416026868585492229216524609739211131993326633970334388991812593549702868877534701822990946125111761892723042376117665640296993581745994557803052315791392349639065203872505E-90";
	coeffs_120[115] = "1.18525924288117432386770939895670573772658621857195305986011196724304231598127227408839423385042572374412446842112646168302015480830234100570192462192015131968307084609177540911503689228342834030959242458698413980031135644018348590823980902427540799814E-91";
	coeffs_120[116] = "-8.5714961216566153236700116412888006837408819915951896129362859520462766617634320531162919426026429378433105901035364956643086394331335747930198070611009941831387116980941022864465946989065467218665543814574849964435089931072761832853235509961870476035E-93";
	coeffs_120[117] = "4.5681983751743456413033268196376305093509590040595182930261094908859252761697530924655649930852283295534503341542929581967081012867692190108698698006237799801339418962091877730207560007839789937153876806052229193448161273005984514504886230869730232561E-94";
	coeffs_120[118] = "-1.5943139155457706045530478744891549581317663177038648406493256399589001327414318955746453934207742828511041930090849236963271943244329753764497401819704943705370596846318480510254313447057477914171472190541408193443142906466279172123681623644325254209E-95";
	coeffs_120[119] = "2.7319125666863032595604997603472305262880292377469053594326527505796348018540179196191192420176181194669607935656210005192217186286873953583571180312679155204061051208771126804209623533044988888808754656646355388901404252058383561064953226611421609762E-97";
	coeffs[3].swap(coeffs_120);
}


/** The Gamma function.
 *  Use the Lanczos approximation. If the coefficients used here are not
 *  sufficiently many or sufficiently accurate, more can be calculated
 *  using the program doc/examples/lanczos.cpp. In that case, be sure to
 *  read the comments in that file. */
const numeric lgamma(const numeric &x)
{
	lanczos_coeffs lc;
	if (lc.sufficiently_accurate(Digits)) {
		cln::cl_N pi_val = cln::pi(cln::default_float_format);
		if (x.real() < 0.5)
			return log(pi_val) - log(sin(pi_val*x.to_cl_N()))
				- lgamma(numeric(1).sub(x));
		cln::cl_N A = lc.calc_lanczos_A(x.to_cl_N());
		cln::cl_N temp = x.to_cl_N() + lc.get_order() - cln::cl_N(1)/2;
   	cln::cl_N result = log(cln::cl_I(2)*pi_val)/2
		              + (x.to_cl_N()-cln::cl_N(1)/2)*log(temp)
		              - temp
		              + log(A);
   	return result;
	}
	else 
		throw dunno();
}

const numeric tgamma(const numeric &x)
{
	lanczos_coeffs lc;
	if (lc.sufficiently_accurate(Digits)) {
		cln::cl_N pi_val = cln::pi(cln::default_float_format);
		if (x.real() < 0.5)
			return pi_val/(sin(pi_val*x))/(tgamma(numeric(1).sub(x)).to_cl_N());
		cln::cl_N A = lc.calc_lanczos_A(x.to_cl_N());
		cln::cl_N temp = x.to_cl_N() + lc.get_order() - cln::cl_N(1)/2;
   	cln::cl_N result
			= sqrt(cln::cl_I(2)*pi_val) * expt(temp, x.to_cl_N()-cln::cl_N(1)/2)
			  * exp(-temp) * A;
   	return result;
	}
	else
		throw dunno();
}


/** The psi function (aka polygamma function).
 *  This is only a stub! */
const numeric psi(const numeric &x)
{
	throw dunno();
}


/** The psi functions (aka polygamma functions).
 *  This is only a stub! */
const numeric psi(const numeric &n, const numeric &x)
{
	throw dunno();
}


/** Factorial combinatorial function.
 *
 *  @param n  integer argument >= 0
 *  @exception range_error (argument must be integer >= 0) */
const numeric factorial(const numeric &n)
{
	if (!n.is_nonneg_integer())
		throw std::range_error("numeric::factorial(): argument must be integer >= 0");
	return numeric(cln::factorial(n.to_int()));
}


/** The double factorial combinatorial function.  (Scarcely used, but still
 *  useful in cases, like for exact results of tgamma(n+1/2) for instance.)
 *
 *  @param n  integer argument >= -1
 *  @return n!! == n * (n-2) * (n-4) * ... * ({1|2}) with 0!! == (-1)!! == 1
 *  @exception range_error (argument must be integer >= -1) */
const numeric doublefactorial(const numeric &n)
{
	if (n.is_equal(*_num_1_p))
		return *_num1_p;
	
	if (!n.is_nonneg_integer())
		throw std::range_error("numeric::doublefactorial(): argument must be integer >= -1");
	
	return numeric(cln::doublefactorial(n.to_int()));
}


/** The Binomial coefficients.  It computes the binomial coefficients.  For
 *  integer n and k and positive n this is the number of ways of choosing k
 *  objects from n distinct objects.  If n is negative, the formula
 *  binomial(n,k) == (-1)^k*binomial(k-n-1,k) is used to compute the result. */
const numeric binomial(const numeric &n, const numeric &k)
{
	if (n.is_integer() && k.is_integer()) {
		if (n.is_nonneg_integer()) {
			if (k.compare(n)!=1 && k.compare(*_num0_p)!=-1)
				return numeric(cln::binomial(n.to_int(),k.to_int()));
			else
				return *_num0_p;
		} else {
			return _num_1_p->power(k)*binomial(k-n-(*_num1_p),k);
		}
	}
	
	// should really be gamma(n+1)/gamma(k+1)/gamma(n-k+1) or a suitable limit
	throw std::range_error("numeric::binomial(): don't know how to evaluate that.");
}


/** Bernoulli number.  The nth Bernoulli number is the coefficient of x^n/n!
 *  in the expansion of the function x/(e^x-1).
 *
 *  @return the nth Bernoulli number (a rational number).
 *  @exception range_error (argument must be integer >= 0) */
const numeric bernoulli(const numeric &nn)
{
	if (!nn.is_integer() || nn.is_negative())
		throw std::range_error("numeric::bernoulli(): argument must be integer >= 0");

	// Method:
	//
	// The Bernoulli numbers are rational numbers that may be computed using
	// the relation
	//
	//     B_n = - 1/(n+1) * sum_{k=0}^{n-1}(binomial(n+1,k)*B_k)
	//
	// with B(0) = 1.  Since the n'th Bernoulli number depends on all the
	// previous ones, the computation is necessarily very expensive.  There are
	// several other ways of computing them, a particularly good one being
	// cl_I s = 1;
	// cl_I c = n+1;
	// cl_RA Bern = 0;
	// for (unsigned i=0; i<n; i++) {
	//     c = exquo(c*(i-n),(i+2));
	//     Bern = Bern + c*s/(i+2);
	//     s = s + expt_pos(cl_I(i+2),n);
	// }
	// return Bern;
	// 
	// But if somebody works with the n'th Bernoulli number she is likely to
	// also need all previous Bernoulli numbers. So we need a complete remember
	// table and above divide and conquer algorithm is not suited to build one
	// up.  The formula below accomplishes this.  It is a modification of the
	// defining formula above but the computation of the binomial coefficients
	// is carried along in an inline fashion.  It also honors the fact that
	// B_n is zero when n is odd and greater than 1.
	// 
	// (There is an interesting relation with the tangent polynomials described
	// in `Concrete Mathematics', which leads to a program a little faster as
	// our implementation below, but it requires storing one such polynomial in
	// addition to the remember table.  This doubles the memory footprint so
	// we don't use it.)

	const unsigned n = nn.to_int();

	// the special cases not covered by the algorithm below
	if (n & 1)
		return (n==1) ? (*_num_1_2_p) : (*_num0_p);
	if (!n)
		return *_num1_p;

	// store nonvanishing Bernoulli numbers here
	static std::vector< cln::cl_RA > results;
	static unsigned next_r = 0;

	// algorithm not applicable to B(2), so just store it
	if (!next_r) {
		results.push_back(cln::recip(cln::cl_RA(6)));
		next_r = 4;
	}
	if (n<next_r)
		return results[n/2-1];

	results.reserve(n/2);
	for (unsigned p=next_r; p<=n;  p+=2) {
		cln::cl_I  c = 1;  // seed for binonmial coefficients
		cln::cl_RA b = cln::cl_RA(p-1)/-2;
		// The CLN manual says: "The conversion from `unsigned int' works only
		// if the argument is < 2^29" (This is for 32 Bit machines. More
		// generally, cl_value_len is the limiting exponent of 2. We must make
		// sure that no intermediates are created which exceed this value. The
		// largest intermediate is (p+3-2*k)*(p/2-k+1) <= (p^2+p)/2.
		if (p < (1UL<<cl_value_len/2)) {
			for (unsigned k=1; k<=p/2-1; ++k) {
				c = cln::exquo(c * ((p+3-2*k) * (p/2-k+1)), (2*k-1)*k);
				b = b + c*results[k-1];
			}
		} else {
			for (unsigned k=1; k<=p/2-1; ++k) {
				c = cln::exquo((c * (p+3-2*k)) * (p/2-k+1), cln::cl_I(2*k-1)*k);
				b = b + c*results[k-1];
			}
 		}
		results.push_back(-b/(p+1));
	}
	next_r = n+2;
	return results[n/2-1];
}


/** Fibonacci number.  The nth Fibonacci number F(n) is defined by the
 *  recurrence formula F(n)==F(n-1)+F(n-2) with F(0)==0 and F(1)==1.
 *
 *  @param n an integer
 *  @return the nth Fibonacci number F(n) (an integer number)
 *  @exception range_error (argument must be an integer) */
const numeric fibonacci(const numeric &n)
{
	if (!n.is_integer())
		throw std::range_error("numeric::fibonacci(): argument must be integer");
	// Method:
	//
	// The following addition formula holds:
	//
	//      F(n+m)   = F(m-1)*F(n) + F(m)*F(n+1)  for m >= 1, n >= 0.
	//
	// (Proof: For fixed m, the LHS and the RHS satisfy the same recurrence
	// w.r.t. n, and the initial values (n=0, n=1) agree. Hence all values
	// agree.)
	// Replace m by m+1:
	//      F(n+m+1) = F(m)*F(n) + F(m+1)*F(n+1)      for m >= 0, n >= 0
	// Now put in m = n, to get
	//      F(2n) = (F(n+1)-F(n))*F(n) + F(n)*F(n+1) = F(n)*(2*F(n+1) - F(n))
	//      F(2n+1) = F(n)^2 + F(n+1)^2
	// hence
	//      F(2n+2) = F(n+1)*(2*F(n) + F(n+1))
	if (n.is_zero())
		return *_num0_p;
	if (n.is_negative())
		if (n.is_even())
			return -fibonacci(-n);
		else
			return fibonacci(-n);
	
	cln::cl_I u(0);
	cln::cl_I v(1);
	cln::cl_I m = cln::the<cln::cl_I>(n.to_cl_N()) >> 1L;  // floor(n/2);
	for (uintL bit=cln::integer_length(m); bit>0; --bit) {
		// Since a squaring is cheaper than a multiplication, better use
		// three squarings instead of one multiplication and two squarings.
		cln::cl_I u2 = cln::square(u);
		cln::cl_I v2 = cln::square(v);
		if (cln::logbitp(bit-1, m)) {
			v = cln::square(u + v) - u2;
			u = u2 + v2;
		} else {
			u = v2 - cln::square(v - u);
			v = u2 + v2;
		}
	}
	if (n.is_even())
		// Here we don't use the squaring formula because one multiplication
		// is cheaper than two squarings.
		return u * ((v << 1) - u);
	else
		return cln::square(u) + cln::square(v);    
}


/** Absolute value. */
const numeric abs(const numeric& x)
{
	return cln::abs(x.to_cl_N());
}


/** Modulus (in positive representation).
 *  In general, mod(a,b) has the sign of b or is zero, and rem(a,b) has the
 *  sign of a or is zero. This is different from Maple's modp, where the sign
 *  of b is ignored. It is in agreement with Mathematica's Mod.
 *
 *  @return a mod b in the range [0,abs(b)-1] with sign of b if both are
 *  integer, 0 otherwise. */
const numeric mod(const numeric &a, const numeric &b)
{
	if (a.is_integer() && b.is_integer())
		return cln::mod(cln::the<cln::cl_I>(a.to_cl_N()),
		                cln::the<cln::cl_I>(b.to_cl_N()));
	else
		return *_num0_p;
}


/** Modulus (in symmetric representation).
 *  Equivalent to Maple's mods.
 *
 *  @return a mod b in the range [-iquo(abs(b)-1,2), iquo(abs(b),2)]. */
const numeric smod(const numeric &a, const numeric &b)
{
	if (a.is_integer() && b.is_integer()) {
		const cln::cl_I b2 = cln::ceiling1(cln::the<cln::cl_I>(b.to_cl_N()) >> 1) - 1;
		return cln::mod(cln::the<cln::cl_I>(a.to_cl_N()) + b2,
		                cln::the<cln::cl_I>(b.to_cl_N())) - b2;
	} else
		return *_num0_p;
}


/** Numeric integer remainder.
 *  Equivalent to Maple's irem(a,b) as far as sign conventions are concerned.
 *  In general, mod(a,b) has the sign of b or is zero, and irem(a,b) has the
 *  sign of a or is zero.
 *
 *  @return remainder of a/b if both are integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric irem(const numeric &a, const numeric &b)
{
	if (b.is_zero())
		throw std::overflow_error("numeric::irem(): division by zero");
	if (a.is_integer() && b.is_integer())
		return cln::rem(cln::the<cln::cl_I>(a.to_cl_N()),
		                cln::the<cln::cl_I>(b.to_cl_N()));
	else
		return *_num0_p;
}


/** Numeric integer remainder.
 *  Equivalent to Maple's irem(a,b,'q') it obeyes the relation
 *  irem(a,b,q) == a - q*b.  In general, mod(a,b) has the sign of b or is zero,
 *  and irem(a,b) has the sign of a or is zero.
 *
 *  @return remainder of a/b and quotient stored in q if both are integer,
 *  0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric irem(const numeric &a, const numeric &b, numeric &q)
{
	if (b.is_zero())
		throw std::overflow_error("numeric::irem(): division by zero");
	if (a.is_integer() && b.is_integer()) {
		const cln::cl_I_div_t rem_quo = cln::truncate2(cln::the<cln::cl_I>(a.to_cl_N()),
		                                               cln::the<cln::cl_I>(b.to_cl_N()));
		q = rem_quo.quotient;
		return rem_quo.remainder;
	} else {
		q = *_num0_p;
		return *_num0_p;
	}
}


/** Numeric integer quotient.
 *  Equivalent to Maple's iquo as far as sign conventions are concerned.
 *  
 *  @return truncated quotient of a/b if both are integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric iquo(const numeric &a, const numeric &b)
{
	if (b.is_zero())
		throw std::overflow_error("numeric::iquo(): division by zero");
	if (a.is_integer() && b.is_integer())
		return cln::truncate1(cln::the<cln::cl_I>(a.to_cl_N()),
	                          cln::the<cln::cl_I>(b.to_cl_N()));
	else
		return *_num0_p;
}


/** Numeric integer quotient.
 *  Equivalent to Maple's iquo(a,b,'r') it obeyes the relation
 *  r == a - iquo(a,b,r)*b.
 *
 *  @return truncated quotient of a/b and remainder stored in r if both are
 *  integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric iquo(const numeric &a, const numeric &b, numeric &r)
{
	if (b.is_zero())
		throw std::overflow_error("numeric::iquo(): division by zero");
	if (a.is_integer() && b.is_integer()) {
		const cln::cl_I_div_t rem_quo = cln::truncate2(cln::the<cln::cl_I>(a.to_cl_N()),
		                                               cln::the<cln::cl_I>(b.to_cl_N()));
		r = rem_quo.remainder;
		return rem_quo.quotient;
	} else {
		r = *_num0_p;
		return *_num0_p;
	}
}


/** Greatest Common Divisor.
 *   
 *  @return  The GCD of two numbers if both are integer, a numerical 1
 *  if they are not. */
const numeric gcd(const numeric &a, const numeric &b)
{
	if (a.is_integer() && b.is_integer())
		return cln::gcd(cln::the<cln::cl_I>(a.to_cl_N()),
		                cln::the<cln::cl_I>(b.to_cl_N()));
	else
		return *_num1_p;
}


/** Least Common Multiple.
 *   
 *  @return  The LCM of two numbers if both are integer, the product of those
 *  two numbers if they are not. */
const numeric lcm(const numeric &a, const numeric &b)
{
	if (a.is_integer() && b.is_integer())
		return cln::lcm(cln::the<cln::cl_I>(a.to_cl_N()),
		                cln::the<cln::cl_I>(b.to_cl_N()));
	else
		return a.mul(b);
}


/** Numeric square root.
 *  If possible, sqrt(x) should respect squares of exact numbers, i.e. sqrt(4)
 *  should return integer 2.
 *
 *  @param x numeric argument
 *  @return square root of x. Branch cut along negative real axis, the negative
 *  real axis itself where imag(x)==0 and real(x)<0 belongs to the upper part
 *  where imag(x)>0. */
const numeric sqrt(const numeric &x)
{
	return cln::sqrt(x.to_cl_N());
}


/** Integer numeric square root. */
const numeric isqrt(const numeric &x)
{
	if (x.is_integer()) {
		cln::cl_I root;
		cln::isqrt(cln::the<cln::cl_I>(x.to_cl_N()), &root);
		return root;
	} else
		return *_num0_p;
}


/** Floating point evaluation of Archimedes' constant Pi. */
ex PiEvalf()
{ 
	return numeric(cln::pi(cln::default_float_format));
}


/** Floating point evaluation of Euler's constant gamma. */
ex EulerEvalf()
{ 
	return numeric(cln::eulerconst(cln::default_float_format));
}


/** Floating point evaluation of Catalan's constant. */
ex CatalanEvalf()
{
	return numeric(cln::catalanconst(cln::default_float_format));
}


/** _numeric_digits default ctor, checking for singleton invariance. */
_numeric_digits::_numeric_digits()
  : digits(17)
{
	// It initializes to 17 digits, because in CLN float_format(17) turns out
	// to be 61 (<64) while float_format(18)=65.  The reason is we want to
	// have a cl_LF instead of cl_SF, cl_FF or cl_DF.
	if (too_late)
		throw(std::runtime_error("I told you not to do instantiate me!"));
	too_late = true;
	cln::default_float_format = cln::float_format(17);

	// add callbacks for built-in functions
	// like ... add_callback(Li_lookuptable);
}


/** Assign a native long to global Digits object. */
_numeric_digits& _numeric_digits::operator=(long prec)
{
	long digitsdiff = prec - digits;
	digits = prec;
	cln::default_float_format = cln::float_format(prec);

	// call registered callbacks
	std::vector<digits_changed_callback>::const_iterator it = callbacklist.begin(),	end = callbacklist.end();
	for (; it != end; ++it) {
		(*it)(digitsdiff);
	}

	return *this;
}


/** Convert global Digits object to native type long. */
_numeric_digits::operator long()
{
	// BTW, this is approx. unsigned(cln::default_float_format*0.301)-1
	return (long)digits;
}


/** Append global Digits object to ostream. */
void _numeric_digits::print(std::ostream &os) const
{
	os << digits;
}


/** Add a new callback function. */
void _numeric_digits::add_callback(digits_changed_callback callback)
{
	callbacklist.push_back(callback);
}


std::ostream& operator<<(std::ostream &os, const _numeric_digits &e)
{
	e.print(os);
	return os;
}

//////////
// static member variables
//////////

// private

bool _numeric_digits::too_late = false;


/** Accuracy in decimal digits.  Only object of this type!  Can be set using
 *  assignment from C++ unsigned ints and evaluated like any built-in type. */
_numeric_digits Digits;

} // namespace GiNaC
