/** @file numeric.h
 *
 *  Makes the interface to the underlying bignum package available. */

/*
 *  This is a modified version of code included with Ginac.  
 *  The modifications and modified version is:
 * 
 *      GiNaC-SAGE Copyright (C) 2008 William Stein
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


/*  The original copyright:
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

#ifndef __GINAC_NUMERIC_H__
#define __GINAC_NUMERIC_H__

#include "basic.h"
#include "constant.h"
#include "ex.h"
#include "py_funcs.h"

#include <gmp.h>
#include <stdexcept>
#include <vector>
#include <iostream>

void ginac_pyinit_Integer(PyObject*);
void ginac_pyinit_Float(PyObject*);
void ginac_pyinit_I(PyObject*);

namespace GiNaC {

enum Type {
//	LONG,
	DOUBLE=1,
	PYOBJECT,
	MPZ,
	MPQ,
//	MPFR,
//	MPFC,
//	MPQC
//	LONG,
};

union Value {
	double _double;
	mpz_t _bigint;
	mpq_t _bigrat;
	PyObject* _pyobject;
};

/** Exception class thrown when a singularity is encountered. */
class pole_error : public std::domain_error {
public:
	explicit pole_error(const std::string& what_arg, int degree);
	int degree() const;
private:
	int deg;
};

/** This class is a wrapper around CLN-numbers within the GiNaC class
 *  hierarchy. Objects of this type may directly be created by the user.*/
class numeric : public basic {
	GINAC_DECLARE_REGISTERED_CLASS(numeric, basic)
	// member functions

	// other constructors
public:
	numeric(const numeric&);
	numeric(int i);
	numeric(unsigned int i);
	numeric(long i);
	numeric(unsigned long i);
	numeric(long numer, long denom);
	numeric(double d);
	numeric(mpz_t bigint);
	numeric(mpq_t bigrat);
	numeric(PyObject*, bool=false);

	~numeric();

	friend std::ostream& operator<<(std::ostream& os, const numeric& s);
	friend void coerce(numeric& new_left, numeric& new_right, const numeric& left, const numeric& right);
	// functions overriding virtual functions from base classes
public:

	int compare_same_type(const numeric& right) const;
	unsigned precedence() const
	{
		return 30;
	}
	bool info(unsigned inf) const;
	bool is_polynomial(const ex & var) const;
	int degree(const ex & s) const;
	int ldegree(const ex & s) const;
	ex coeff(const ex & s, int n = 1) const;
	bool has(const ex &other, unsigned options = 0) const;
	ex eval(int level = 0) const;
	ex evalf(int level = 0, PyObject* parent = nullptr) const;

	ex subs(const exmap & m, unsigned options = 0) const
	{
		return subs_one_level(m, options);
	} // overwrites basic::subs() for performance reasons
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0) const;
	ex to_rational(exmap & repl) const;
	ex to_polynomial(exmap & repl) const;
	numeric integer_content() const;
	numeric max_coefficient() const;
	ex conjugate() const;
	ex real_part() const;
	ex imag_part() const;
protected:

	/** Implementation of ex::diff for a numeric always returns 0.
	 *  @see ex::diff */
	ex derivative(const symbol &s) const
	{
		return 0;
	}
	bool is_equal_same_type(const basic &other) const;
	int64_t calchash() const;

	// new virtual functions which can be overridden by derived classes
	// (none)

	// non-virtual functions in this class
public:
	const numeric add(const numeric &other) const;
	const numeric sub(const numeric &other) const;
	const numeric mul(const numeric &other) const;
	const numeric div(const numeric &other) const;
	const numeric power(const numeric &other) const;
	const numeric & add_dyn(const numeric &other) const;
	const numeric & sub_dyn(const numeric &other) const;
	const numeric & mul_dyn(const numeric &other) const;
	const numeric & div_dyn(const numeric &other) const;
	const numeric & power_dyn(const numeric &other) const;
	const numeric & operator=(int i);
	const numeric & operator=(unsigned int i);
	const numeric & operator=(long i);
	const numeric & operator=(unsigned long i);
	const numeric & operator=(double d);
	const numeric & operator=(const numeric& x);
	const numeric negative() const;
	const numeric inverse() const;
	const numeric step() const;
	int csgn() const;
	int compare(const numeric &other) const;
	bool is_equal(const numeric &other) const;
	bool is_zero() const;
	bool is_positive() const;
	bool is_negative() const;
	bool is_integer() const;
	bool is_pos_integer() const;
	bool is_nonneg_integer() const;
	bool is_even() const;
	bool is_odd() const;
	bool is_prime() const;
	bool is_rational() const;
	bool is_real() const;
	bool is_cinteger() const;
	bool is_crational() const;
	bool is_exact() const;
	bool operator==(const numeric &other) const;
	bool operator!=(const numeric &other) const;
	bool operator<(const numeric &other) const;
	bool operator<=(const numeric &other) const;
	bool operator>(const numeric &other) const;
	bool operator>=(const numeric &other) const;
	bool is_parent_pos_char() const;
	int get_parent_char() const;
	int to_int() const
	{
		return (int)to_long();
	}
	long to_long() const;
	double to_double() const;
	PyObject* to_pyobject() const;
	const numeric real() const;
	const numeric imag() const;
	const numeric numer() const;
	const numeric denom() const;

	const numeric exp() const;
	const numeric log() const;
	const numeric sin() const;
	const numeric cos() const;
	const numeric tan() const;
	const numeric asin() const;
	const numeric acos() const;
	const numeric atan() const;
	const numeric atan(const numeric &y) const;
	const numeric sinh() const;
	const numeric cosh() const;
	const numeric tanh() const;
	const numeric asinh() const;
	const numeric acosh() const;
	const numeric atanh() const;
	const numeric Li2() const;
	const numeric Li2(const numeric &n, PyObject* parent) const;
	const numeric zeta() const;
	const numeric lgamma() const;
	const numeric tgamma() const;
	const numeric psi() const;
	const numeric psi(const numeric &n) const;
	const numeric factorial() const;
	const numeric doublefactorial() const;
	const numeric binomial(const numeric &k) const;
	const numeric bernoulli() const;
	const numeric fibonacci() const;
	const numeric isqrt() const;
	const numeric sqrt() const;
	const numeric abs() const;
	const numeric mod(const numeric &b) const;
	const numeric _smod(const numeric &b) const;
	ex smod(const numeric &b) const;
	const numeric irem(const numeric &b) const;
	const numeric iquo(const numeric &b) const;
	const numeric iquo(const numeric &b, numeric &r) const;
	const numeric gcd(const numeric &b) const;
	const numeric lcm(const numeric &b) const;
	
	int int_length() const;

protected:
	void print_numeric(const print_context & c, const char *par_open,
		const char *par_close, const char *imag_sym,
		const char *mul_sym, unsigned level, bool latex) const;
	void do_print(const print_context & c, unsigned level) const;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_csrc(const print_csrc & c, unsigned level) const;
	void do_print_tree(const print_tree & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const;

//	numeric operator()(const int& x);

	// member variables

protected:
    Type t;
    Value v;
    int64_t hash;
    bool is_hashable = true;
};


// global constants

extern numeric I;
extern PyObject *RR, *CC;

// global functions

void coerce(numeric& new_left, numeric& new_right, const numeric& left, const numeric& right);

const numeric exp(const numeric &x);
const numeric log(const numeric &x);
const numeric sin(const numeric &x);
const numeric cos(const numeric &x);
const numeric tan(const numeric &x);
const numeric asin(const numeric &x);
const numeric acos(const numeric &x);
const numeric atan(const numeric &x);
const numeric atan(const numeric &y, const numeric &x);
const numeric sinh(const numeric &x);
const numeric cosh(const numeric &x);
const numeric tanh(const numeric &x);
const numeric asinh(const numeric &x);
const numeric acosh(const numeric &x);
const numeric atanh(const numeric &x);
const numeric Li2(const numeric &x);
const numeric Li2(const numeric &x, const numeric &n, PyObject* parent);
const numeric zeta(const numeric &x);
const numeric lgamma(const numeric &x);
const numeric tgamma(const numeric &x);
const numeric psi(const numeric &x);
const numeric psi(const numeric &n, const numeric &x);
const numeric factorial(const numeric &n);
const numeric doublefactorial(const numeric &n);
const numeric binomial(const numeric &n, const numeric &k);
const numeric bernoulli(const numeric &n);
const numeric fibonacci(const numeric &n);
const numeric isqrt(const numeric &x);
const numeric sqrt(const numeric &x);
const numeric abs(const numeric &x);
const numeric mod(const numeric &a, const numeric &b);
const numeric smod(const numeric &a, const numeric &b);
const numeric irem(const numeric &a, const numeric &b);
const numeric iquo(const numeric &a, const numeric &b);
const numeric iquo(const numeric &a, const numeric &b, numeric &r);
const numeric gcd(const numeric &a, const numeric &b);
const numeric lcm(const numeric &a, const numeric &b);

// wrapper functions around member functions

inline const numeric pow(const numeric &x, const numeric &y)
{
	return x.power(y);
}

inline const numeric inverse(const numeric &x)
{
	return x.inverse();
}

inline numeric step(const numeric &x)
{
	return x.step();
}

inline int csgn(const numeric &x)
{
	return x.csgn();
}

inline bool is_zero(const numeric &x)
{
	return x.is_zero();
}

inline bool is_positive(const numeric &x)
{
	return x.is_positive();
}

inline bool is_negative(const numeric &x)
{
	return x.is_negative();
}

inline bool is_integer(const numeric &x)
{
	return x.is_integer();
}

inline bool is_pos_integer(const numeric &x)
{
	return x.is_pos_integer();
}

inline bool is_nonneg_integer(const numeric &x)
{
	return x.is_nonneg_integer();
}

inline bool is_even(const numeric &x)
{
	return x.is_even();
}

inline bool is_odd(const numeric &x)
{
	return x.is_odd();
}

inline bool is_prime(const numeric &x)
{
	return x.is_prime();
}

inline bool is_rational(const numeric &x)
{
	return x.is_rational();
}

inline bool is_real(const numeric &x)
{
	return x.is_real();
}

inline bool is_cinteger(const numeric &x)
{
	return x.is_cinteger();
}

inline bool is_crational(const numeric &x)
{
	return x.is_crational();
}

inline int to_int(const numeric &x)
{
	return x.to_int();
}

inline long to_long(const numeric &x)
{
	return x.to_long();
}

inline double to_double(const numeric &x)
{
	return x.to_double();
}

inline const numeric real(const numeric &x)
{
	return x.real();
}

inline const numeric imag(const numeric &x)
{
	return x.imag();
}

inline const numeric numer(const numeric &x)
{
	return x.numer();
}

inline const numeric denom(const numeric &x)
{
	return x.denom();
}

// numeric evaluation functions for class constant objects:

ex ConstantEvalf(unsigned serial, PyObject* parent = nullptr);
ex UnsignedInfinityEvalf(unsigned serial, PyObject* parent = nullptr);
ex InfinityEvalf(unsigned serial, PyObject* parent = nullptr);
ex NegInfinityEvalf(unsigned serial, PyObject* parent = nullptr);

} // namespace GiNaC

#endif // ndef __GINAC_NUMERIC_H__
