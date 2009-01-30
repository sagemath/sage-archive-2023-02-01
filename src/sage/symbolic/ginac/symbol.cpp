/** @file symbol.cpp
 *
 *  Implementation of GiNaC's symbolic objects. */

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

#include <string>
#include <stdexcept>

#include "symbol.h"
#include "power.h"
#include "mul.h"
#include "lst.h"
#include "archive.h"
#include "tostring.h"
#include "utils.h"
#include "inifcns.h"

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(symbol, basic,
  print_func<print_context>(&symbol::do_print).
  print_func<print_latex>(&symbol::do_print_latex).
  print_func<print_tree>(&symbol::do_print_tree).
  print_func<print_python_repr>(&symbol::do_print_python_repr))

//////////
// default constructor
//////////

// symbol

symbol::symbol()
 : inherited(&symbol::tinfo_static), asexinfop(new assigned_ex_info), serial(next_serial++), name(autoname_prefix() + ToString(serial)), TeX_name(name), domain(domain::complex), ret_type(return_types::commutative), ret_type_tinfo(&symbol::tinfo_static)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

// realsymbol

realsymbol::realsymbol()
{
	domain = domain::real;
}

// possymbol

possymbol::possymbol()
{
	domain = domain::positive;
}

//////////
// other constructors
//////////

// public

// symbol

symbol::symbol(const std::string & initname, unsigned domain)
 : inherited(&symbol::tinfo_static), asexinfop(new assigned_ex_info), serial(next_serial++), name(initname), TeX_name(default_TeX_name()), domain(domain), ret_type(return_types::commutative), ret_type_tinfo(&symbol::tinfo_static)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

symbol::symbol(const std::string & initname, unsigned rt, tinfo_t rtt, unsigned domain)
 : inherited(&symbol::tinfo_static), asexinfop(new assigned_ex_info), serial(next_serial++), name(initname), TeX_name(default_TeX_name()), domain(domain), ret_type(rt), ret_type_tinfo(rtt)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

symbol::symbol(const std::string & initname, const std::string & texname, unsigned domain)
 : inherited(&symbol::tinfo_static), asexinfop(new assigned_ex_info), serial(next_serial++), name(initname), TeX_name(texname), domain(domain), ret_type(return_types::commutative), ret_type_tinfo(&symbol::tinfo_static)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

symbol::symbol(const std::string & initname, const std::string & texname, unsigned rt, tinfo_t rtt, unsigned domain)
 : inherited(&symbol::tinfo_static), asexinfop(new assigned_ex_info), serial(next_serial++), name(initname), TeX_name(texname), domain(domain), ret_type(rt), ret_type_tinfo(rtt)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

// realsymbol
	
realsymbol::realsymbol(const std::string & initname, unsigned domain)
 : symbol(initname, domain) { }

realsymbol::realsymbol(const std::string & initname, const std::string & texname, unsigned domain)
 : symbol(initname, texname, domain) { }

realsymbol::realsymbol(const std::string & initname, unsigned rt, tinfo_t rtt, unsigned domain)
 : symbol(initname, rt, rtt, domain) { }

realsymbol::realsymbol(const std::string & initname, const std::string & texname, unsigned rt, tinfo_t rtt, unsigned domain)
 : symbol(initname, texname, rt, rtt, domain) { }

// possymbol
	
possymbol::possymbol(const std::string & initname, unsigned domain)
 : symbol(initname, domain) { }

possymbol::possymbol(const std::string & initname, const std::string & texname, unsigned domain)
 : symbol(initname, texname, domain) { }

possymbol::possymbol(const std::string & initname, unsigned rt, tinfo_t rtt, unsigned domain)
 : symbol(initname, rt, rtt, domain) { }

possymbol::possymbol(const std::string & initname, const std::string & texname, unsigned rt, tinfo_t rtt, unsigned domain)
 : symbol(initname, texname, rt, rtt, domain) { }

//////////
// archiving
//////////

/** Construct object from archive_node. */
symbol::symbol(const archive_node &n, lst &sym_lst)
 : inherited(n, sym_lst), asexinfop(new assigned_ex_info), serial(next_serial++)
{
	if (!n.find_string("name", name))
		name = autoname_prefix() + ToString(serial);
	if (!n.find_string("TeXname", TeX_name))
		TeX_name = default_TeX_name();
	if (!n.find_unsigned("domain", domain))
		domain = domain::complex;
	if (!n.find_unsigned("return_type", ret_type))
		ret_type = return_types::commutative;
	setflag(status_flags::evaluated | status_flags::expanded);
}

/** Unarchive the object. */
ex symbol::unarchive(const archive_node &n, lst &sym_lst)
{
	ex s = (new symbol(n, sym_lst))->setflag(status_flags::dynallocated);

	// If symbol is in sym_lst, return the existing symbol
	for (lst::const_iterator it = sym_lst.begin(); it != sym_lst.end(); ++it) {
		if (is_a<symbol>(*it) && (ex_to<symbol>(*it).name == ex_to<symbol>(s).name))
			return *it;
	}

	// Otherwise add new symbol to list and return it
	sym_lst.append(s);
	return s;
}

/** Archive the object. */
void symbol::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_string("name", name);
	if (TeX_name != default_TeX_name())
		n.add_string("TeX_name", TeX_name);
	if (domain != domain::complex)
		n.add_unsigned("domain", domain);
	if (ret_type != return_types::commutative)
		n.add_unsigned("return_type", ret_type);
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void symbol::do_print(const print_context & c, unsigned level) const
{
	c.s << name;
}

void symbol::do_print_latex(const print_latex & c, unsigned level) const
{
	c.s << TeX_name;
}

void symbol::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << name << " (" << class_name() << ")" << " @" << this
	    << ", serial=" << serial
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", domain=" << domain
	    << std::endl;
}

void symbol::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << "('" << name;
	if (TeX_name != default_TeX_name())
		c.s << "','" << TeX_name;
	c.s << "')";
}

bool symbol::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::symbol:
		case info_flags::polynomial:
		case info_flags::integer_polynomial: 
		case info_flags::cinteger_polynomial: 
		case info_flags::rational_polynomial: 
		case info_flags::crational_polynomial: 
		case info_flags::rational_function: 
		case info_flags::expanded:
			return true;
		case info_flags::real:
			return domain == domain::real || domain == domain::positive;
		case info_flags::positive:
		case info_flags::nonnegative:
			return domain == domain::positive;
		case info_flags::has_indices:
			return false;
	}
	return inherited::info(inf);
}

ex symbol::eval(int level) const
{
	if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	if (asexinfop->is_assigned) {
		setflag(status_flags::evaluated);
		if (level==1)
			return (asexinfop->assigned_expression);
		else
			return (asexinfop->assigned_expression).eval(level);
	} else {
		return this->hold();
	}
}

ex symbol::conjugate() const
{
	if (this->domain == domain::complex) {
		return conjugate_function(*this).hold();
	} else {
		return *this;
	}
}

ex symbol::real_part() const
{
	if (domain==domain::real || domain==domain::positive)
		return *this;
	return real_part_function(*this).hold();
}

ex symbol::imag_part() const
{
	if (domain==domain::real || domain==domain::positive)
		return 0;
	return imag_part_function(*this).hold();
}

bool symbol::is_polynomial(const ex & var) const
{
	return true;
}

// protected

/** Implementation of ex::diff() for single differentiation of a symbol.
 *  It returns 1 or 0.
 *
 *  @see ex::diff */
ex symbol::derivative(const symbol & s) const
{
	if (compare_same_type(s))
		return _ex0;
	else
		return _ex1;
}

int symbol::compare(const basic& other) const
{
	static const tinfo_t pow_id = find_tinfo_key("power");
	static const tinfo_t mul_id = find_tinfo_key("mul");
	static const tinfo_t function_id = find_tinfo_key("function");
	static const tinfo_t fderivative_id = find_tinfo_key("fderivative");
	const tinfo_t typeid_this = tinfo();
	const tinfo_t typeid_other = other.tinfo();
	if (typeid_this==typeid_other) {
		GINAC_ASSERT(typeid(*this)==typeid(other));
		return compare_same_type(other);
	} else if (typeid_other == pow_id) {
		return -static_cast<const power&>(other).compare_symbol(*this);
	} else if (typeid_other == mul_id) {
		return -static_cast<const mul&>(other).compare_symbol(*this);
	} else if (typeid_other == function_id ||
			typeid_other == fderivative_id) {
		return -1;
	} else {
		return (typeid_this<typeid_other ? -1 : 1);
	}
}

int symbol::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<symbol>(other));
	const symbol *o = static_cast<const symbol *>(&other);
	if (serial==o->serial) return 0;
	
	// SAGE/Pynac: Sorting based on creation order doesn't work for Sage.
	// instead we sort on variable name. -- William Stein

	return name < (o->name) ? -1 : 1;
	
	// This is what Ginac used to return.  It fits with their
	// philosophy that symbols names have no intrinsic meaning,
	// only the order of creation of symbols matters.  Also, they
	// allow multiple symbols with the same name.  Our Pynac
	// wrapper does not allow for this, and we find the behavior
	// way too confusing and ad hoc with this convention. 

	// return serial < o->serial ? -1 : 1;
}

bool symbol::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<symbol>(other));
	const symbol *o = static_cast<const symbol *>(&other);
	return serial==o->serial;
}

unsigned symbol::calchash() const
{
  // SAGE -- instead of some random hash, which results in variables
  // printing and being stored in confusing nondeterministic order as
  // in Ginac, we use a hash designed to make variables print in a
  // sensible order.  In particular we turn the characters of the
  // variable name into an unsigned integer, basically by viewing
  // them as giving a base 256 encoded number. 
  //
  // We do this instead of overloading compare because it is FAST -- the
  // hash is computed once and for all, then all future compares are
  // essentially instant.  This is likely an extremely time-critical part of
  // ginac judging by the benchmarking code in ex.h.
  // 
  // WARNING: We hash strings to 32-bit ints; the result is that
  // variable names up to 4 characters are sorted lexicographically
  // automatically by their hash.  Longer names with same first 4
  // characters result in a hash collision, which of course will take
  // a tiny bit longer.  This is correctly dealt with in the
  // compare_same_type function above.
  // -- William Stein

  hashvalue = 0;
  unsigned int b = 16777216;   // 2^32 / 256
  unsigned int m = (name.length() < 4) ? name.length() : 4;   // = min(4, name.length())
  for (int i=0; i < m; i++) {
    hashvalue += (int)name[i]*b;
    b /= 256;
  }
  setflag(status_flags::hash_calculated);
  return hashvalue;
  
  // Original code
  /*
    hashvalue = golden_ratio_hash((p_int)tinfo() ^ serial);
    setflag(status_flags::hash_calculated);
    return hashvalue;
  */
}

//////////
// virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// public

void symbol::assign(const ex & value)
{
	asexinfop->is_assigned = true;
	asexinfop->assigned_expression = value;
	clearflag(status_flags::evaluated | status_flags::expanded);
}

void symbol::unassign()
{
	if (asexinfop->is_assigned) {
		asexinfop->is_assigned = false;
		asexinfop->assigned_expression = _ex0;
	}
	setflag(status_flags::evaluated | status_flags::expanded);
}

// private

/** Symbols not constructed with a string get one assigned using this
 *  prefix and a number. */
std::string & symbol::autoname_prefix()
{
	static std::string *s = new std::string("symbol");
	return *s;
}

/** Return default TeX name for symbol. This recognizes some greek letters. */
std::string symbol::default_TeX_name() const
{
	if (name=="alpha"        || name=="beta"         || name=="gamma"
	 || name=="delta"        || name=="epsilon"      || name=="varepsilon"
	 || name=="zeta"         || name=="eta"          || name=="theta"
	 || name=="vartheta"     || name=="iota"         || name=="kappa"
	 || name=="lambda"       || name=="mu"           || name=="nu"
	 || name=="xi"           || name=="omicron"      || name=="pi"
	 || name=="varpi"        || name=="rho"          || name=="varrho"
	 || name=="sigma"        || name=="varsigma"     || name=="tau"
	 || name=="upsilon"      || name=="phi"          || name=="varphi"
	 || name=="chi"          || name=="psi"          || name=="omega"
	 || name=="Gamma"        || name=="Delta"        || name=="Theta"
	 || name=="Lambda"       || name=="Xi"           || name=="Pi"
	 || name=="Sigma"        || name=="Upsilon"      || name=="Phi"
	 || name=="Psi"          || name=="Omega")
		return "\\" + name;
	else
		return name;
}

//////////
// static member variables
//////////

// private

unsigned symbol::next_serial = 0;

//////////
// subclass assigned_ex_info
//////////

/** Default ctor.  Defaults to unassigned. */
symbol::assigned_ex_info::assigned_ex_info() throw() : is_assigned(false)
{
}

} // namespace GiNaC
