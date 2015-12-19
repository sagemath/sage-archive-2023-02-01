/** @file constant.cpp
 *
 *  Implementation of GiNaC's constant types and some special constants. */

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

#include "py_funcs.h"
#include "constant.h"
#include "numeric.h"
#include "ex.h"
#include "archive.h"
#include "utils.h"
#include "inifcns.h"

#include <string>
#include <stdexcept>
#include <iostream>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(constant, basic,
  print_func<print_context>(&constant::do_print).
  print_func<print_latex>(&constant::do_print_latex).
  print_func<print_tree>(&constant::do_print_tree).
  print_func<print_python_repr>(&constant::do_print_python_repr))

//////////
// default constructor
//////////

// public

constant::constant() : basic(&constant::tinfo_static), ef(nullptr), serial(next_serial++), domain(domain::complex)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

// public

constant::constant(std::string  initname, evalffunctype efun, const std::string & texname, unsigned dm)
  : basic(&constant::tinfo_static), name(std::move(initname)), ef(efun), serial(next_serial++), domain(dm)
{
	if (texname.empty())
		TeX_name = "\\mbox{" + name + "}";
	else
		TeX_name = texname;
	setflag(status_flags::evaluated | status_flags::expanded);
}

constant::constant(std::string  initname, const numeric & initnumber, const std::string & texname, unsigned dm)
  : basic(&constant::tinfo_static), name(std::move(initname)), ef(nullptr), number(initnumber), serial(next_serial++), domain(dm)
{
	if (texname.empty())
		TeX_name = "\\mbox{" + name + "}";
	else
		TeX_name = texname;
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// archiving
//////////

constant::constant(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst) {}

ex constant::unarchive(const archive_node &n, lst &sym_lst)
{
	// Find constant by name (!! this is bad: 'twould be better if there
	// was a list of all global constants that we could search)
	constant ans;
	std::string s;
	if (n.find_string("name", s)) {
		if (s == Pi.name)
			return Pi;
		else if (s == Catalan.name)
			return Catalan;
		else if (s == Euler.name)
			return Euler;
		else {
			ans = py_funcs.py_get_constant(s.c_str());
			if (PyErr_Occurred()) {
				throw std::runtime_error("error while unarchiving constant");
			}
		}

		
            //throw (std::runtime_error("unknown constant '" + s + "' in archive"));
			return ans;      
	} else
		throw (std::runtime_error("unnamed constant in archive"));
}

void constant::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_string("name", name);
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void constant::do_print(const print_context & c, unsigned level) const
{
	c.s << name;
}

void constant::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << name << " (" << class_name() << ")" << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << std::endl;
}

void constant::do_print_latex(const print_latex & c, unsigned level) const
{
	c.s << TeX_name;
}

void constant::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << "('" << name << "'";
	if (TeX_name != "\\mbox{" + name + "}")
		c.s << ",TeX_name='" << TeX_name << "'";
	c.s << ')';
}

bool constant::info(unsigned inf) const
{
	if (inf == info_flags::polynomial)
		return true;
        if (inf == info_flags::inexact)
                return false;
	if (inf == info_flags::real)
		return domain==domain::real || domain==domain::positive ;
	if (inf==info_flags::positive || inf==info_flags::nonnegative)
		return domain == domain::positive;
	if (inf==info_flags::infinity) {
		return domain == domain::infinity;
	}
	else
		return inherited::info(inf);
}

ex constant::evalf(int level, PyObject* parent) const
{
	if (ef!=nullptr) {
        return ef(serial, parent);
	} else {
		return number.evalf(level, parent);
	}
	return *this;
}

bool constant::is_polynomial(const ex & var) const
{
	return true;
}

ex constant::conjugate() const
{
	if ( domain==domain::real || domain==domain::positive )
		return *this;
	return conjugate_function(*this).hold();
}

ex constant::real_part() const
{
	if ( domain==domain::real || domain==domain::positive )
		return *this;
	return real_part_function(*this).hold();
}

ex constant::imag_part() const
{
	if ( domain==domain::real || domain==domain::positive )
		return 0;
	return imag_part_function(*this).hold();
}

// protected

/** Implementation of ex::diff() for a constant always returns 0.
 *
 *  @see ex::diff */
ex constant::derivative(const symbol & s) const
{
	return _ex0;
}

int constant::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<constant>(other));
	const constant &o = static_cast<const constant &>(other);

	if (serial == o.serial)
		return 0;
	else
		return serial < o.serial ? -1 : 1;
}

bool constant::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<constant>(other));
	const constant &o = static_cast<const constant &>(other);

	return serial == o.serial;
}

int64_t constant::calchash() const
{
	hashvalue = golden_ratio_hash((p_int)tinfo() ^ serial);

	setflag(status_flags::hash_calculated);

	return hashvalue;
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// none

//////////
// static member variables
//////////

unsigned constant::next_serial = 0;

//////////
// global constants
//////////

/**  Pi. (3.14159...) Calls python function py_eval_constant() for evalf(). */
const constant Pi("pi", ConstantEvalf, "\\pi", domain::positive);

/** Euler's constant. (0.57721...)  Sometimes called Euler-Mascheroni constant.
 *  Calls python function py_eval_consant for evalf(). */
const constant Euler("euler_gamma", ConstantEvalf, "\\gamma_E", domain::positive);

/** Catalan's constant. (0.91597...)  
 *  Calls python function py_eval_constant for evalf(). */
const constant Catalan("catalan", ConstantEvalf, "G", domain::positive);

} // namespace GiNaC
