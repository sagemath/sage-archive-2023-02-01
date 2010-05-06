#  This perl script automatically generates function.h and function.cpp

#  function.pl options: \$maxargs=${maxargs}
# 
#  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
# 
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

$maxargs=14;

sub generate_seq {
	my ($seq_template,$n)=@_;
	my ($res,$N);
	
	$res='';
	for ($N=1; $N<=$n; $N++) {
		$res .= eval('"' . $seq_template . '"');
		if ($N!=$n) {
			$res .= ', ';
		}
	}
	return $res;
}

sub generate_from_to {
	my ($template,$seq_template1,$seq_template2,$seq_template3,$from,$to)=@_;
	my ($res,$N,$SEQ);

	$res='';
	for ($N=$from; $N<=$to; $N++) {
		$SEQ1=generate_seq($seq_template1,$N);
		$SEQ2=generate_seq($seq_template2,$N);
		$SEQ3=generate_seq($seq_template3,$N);
		$res .= eval('"' . $template . '"');
		$SEQ1=''; # to avoid main::SEQ1 used only once warning
		$SEQ2=''; # same as above
		$SEQ3=''; # same as above
	}
	return $res;
}

sub generate {
	my ($template,$seq_template1,$seq_template2,$seq_template3)=@_;
	return generate_from_to($template,$seq_template1,$seq_template2,$seq_template3,1,$maxargs);
}

$declare_function_macro = generate(
	<<'END_OF_DECLARE_FUNCTION_MACRO','typename T${N}','const T${N} & p${N}','GiNaC::ex(p${N})');
#define DECLARE_FUNCTION_${N}P(NAME) \\
class NAME##_SERIAL { public: static unsigned serial; }; \\
const unsigned NAME##_NPARAMS = ${N}; \\
template<${SEQ1}> const GiNaC::function NAME(${SEQ2}) { \\
	return GiNaC::function(NAME##_SERIAL::serial, ${SEQ3}); \\
}

END_OF_DECLARE_FUNCTION_MACRO

$typedef_eval_funcp=generate(
'typedef ex (* eval_funcp_${N})(${SEQ1});'."\n",
'const ex &','','');

$typedef_evalf_funcp=generate(
'typedef ex (* evalf_funcp_${N})(${SEQ1});'."\n",
'const ex &','','');

$typedef_conjugate_funcp=generate(
'typedef ex (* conjugate_funcp_${N})(${SEQ1});'."\n",
'const ex &','','');

$typedef_real_part_funcp=generate(
'typedef ex (* real_part_funcp_${N})(${SEQ1});'."\n",
'const ex &','','');

$typedef_imag_part_funcp=generate(
'typedef ex (* imag_part_funcp_${N})(${SEQ1});'."\n",
'const ex &','','');

$typedef_derivative_funcp=generate(
'typedef ex (* derivative_funcp_${N})(${SEQ1}, unsigned);'."\n",
'const ex &','','');

$typedef_power_funcp=generate(
'typedef ex (* power_funcp_${N})(${SEQ1}, const ex &);'."\n",
'const ex &','','');

$typedef_series_funcp=generate(
'typedef ex (* series_funcp_${N})(${SEQ1}, const relational &, int, unsigned);'."\n",
'const ex &','','');

$typedef_print_funcp=generate(
'typedef void (* print_funcp_${N})(${SEQ1}, const print_context &);'."\n",
'const ex &','','');

$eval_func_interface=generate('    function_options & eval_func(eval_funcp_${N} e);'."\n",'','','');

$evalf_func_interface=generate('    function_options & evalf_func(evalf_funcp_${N} ef);'."\n",'','','');

$conjugate_func_interface=generate('    function_options & conjugate_func(conjugate_funcp_${N} d);'."\n",'','','');

$real_part_func_interface=generate('    function_options & real_part_func(real_part_funcp_${N} d);'."\n",'','','');

$imag_part_func_interface=generate('    function_options & imag_part_func(imag_part_funcp_${N} d);'."\n",'','','');

$derivative_func_interface=generate('    function_options & derivative_func(derivative_funcp_${N} d);'."\n",'','','');

$power_func_interface=generate('    function_options & power_func(power_funcp_${N} d);'."\n",'','','');

$series_func_interface=generate('    function_options & series_func(series_funcp_${N} s);'."\n",'','','');

$print_func_interface=generate(
	<<'END_OF_PRINT_FUNC_INTERFACE','','','');
    template <class Ctx> function_options & print_func(print_funcp_${N} p)
    {
    	test_and_set_nparams(${N});
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
END_OF_PRINT_FUNC_INTERFACE

$constructors_interface=generate(
'    function(unsigned ser, ${SEQ1});'."\n",
'const ex & param${N}','','');

$constructors_implementation=generate(
	<<'END_OF_CONSTRUCTORS_IMPLEMENTATION','const ex & param${N}','param${N}','');
function::function(unsigned ser, ${SEQ1})
	: exprseq(${SEQ2}), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
END_OF_CONSTRUCTORS_IMPLEMENTATION

$eval_switch_statement=generate(
	<<'END_OF_EVAL_SWITCH_STATEMENT','seq[${N}-1]','','');
	case ${N}:
		eval_result = ((eval_funcp_${N})(opt.eval_f))(${SEQ1});
		break;
END_OF_EVAL_SWITCH_STATEMENT

$evalf_switch_statement=generate(
	<<'END_OF_EVALF_SWITCH_STATEMENT','eseq[${N}-1]','','');
	case ${N}:
		return ((evalf_funcp_${N})(opt.evalf_f))(${SEQ1});
END_OF_EVALF_SWITCH_STATEMENT

$conjugate_switch_statement=generate(
	<<'END_OF_DIFF_SWITCH_STATEMENT','seq[${N}-1]','','');
	case ${N}:
		return ((conjugate_funcp_${N})(opt.conjugate_f))(${SEQ1});
END_OF_DIFF_SWITCH_STATEMENT

$real_part_switch_statement=generate(
	<<'END_OF_DIFF_SWITCH_STATEMENT','seq[${N}-1]','','');
	case ${N}:
		return ((real_part_funcp_${N})(opt.real_part_f))(${SEQ1});
END_OF_DIFF_SWITCH_STATEMENT

$imag_part_switch_statement=generate(
	<<'END_OF_DIFF_SWITCH_STATEMENT','seq[${N}-1]','','');
	case ${N}:
		return ((imag_part_funcp_${N})(opt.imag_part_f))(${SEQ1});
END_OF_DIFF_SWITCH_STATEMENT

$diff_switch_statement=generate(
	<<'END_OF_DIFF_SWITCH_STATEMENT','seq[${N}-1]','','');
	case ${N}:
		return ((derivative_funcp_${N})(opt.derivative_f))(${SEQ1},diff_param);
END_OF_DIFF_SWITCH_STATEMENT

$power_switch_statement=generate(
	<<'END_OF_POWER_SWITCH_STATEMENT','seq[${N}-1]','','');
	case ${N}:
		return ((power_funcp_${N})(opt.power_f))(${SEQ1},power_param);
END_OF_POWER_SWITCH_STATEMENT

$series_switch_statement=generate(
	<<'END_OF_SERIES_SWITCH_STATEMENT','seq[${N}-1]','','');
	case ${N}:
		try {
			res = ((series_funcp_${N})(opt.series_f))(${SEQ1},r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
END_OF_SERIES_SWITCH_STATEMENT

$print_switch_statement=generate(
	<<'END_OF_PRINT_SWITCH_STATEMENT','seq[${N}-1]','','');
		case ${N}:
			((print_funcp_${N})(pdt[id]))(${SEQ1}, c);
			break;
END_OF_PRINT_SWITCH_STATEMENT

$eval_func_implementation=generate(
	<<'END_OF_EVAL_FUNC_IMPLEMENTATION','','','');
function_options & function_options::eval_func(eval_funcp_${N} e)
{
	test_and_set_nparams(${N});
	eval_f = eval_funcp(e);
	return *this;
}
END_OF_EVAL_FUNC_IMPLEMENTATION

$evalf_func_implementation=generate(
	<<'END_OF_EVALF_FUNC_IMPLEMENTATION','','','');
function_options & function_options::evalf_func(evalf_funcp_${N} ef)
{
	test_and_set_nparams(${N});
	evalf_f = evalf_funcp(ef);
	return *this;
}
END_OF_EVALF_FUNC_IMPLEMENTATION

$conjugate_func_implementation=generate(
	<<'END_OF_CONJUGATE_FUNC_IMPLEMENTATION','','','');
function_options & function_options::conjugate_func(conjugate_funcp_${N} c)
{
	test_and_set_nparams(${N});
	conjugate_f = conjugate_funcp(c);
	return *this;
}
END_OF_CONJUGATE_FUNC_IMPLEMENTATION

$real_part_func_implementation=generate(
	<<'END_OF_REAL_PART_FUNC_IMPLEMENTATION','','','');
function_options & function_options::real_part_func(real_part_funcp_${N} c)
{
	test_and_set_nparams(${N});
	real_part_f = real_part_funcp(c);
	return *this;
}
END_OF_REAL_PART_FUNC_IMPLEMENTATION

$imag_part_func_implementation=generate(
	<<'END_OF_IMAG_PART_FUNC_IMPLEMENTATION','','','');
function_options & function_options::imag_part_func(imag_part_funcp_${N} c)
{
	test_and_set_nparams(${N});
	imag_part_f = imag_part_funcp(c);
	return *this;
}
END_OF_IMAG_PART_FUNC_IMPLEMENTATION

$derivative_func_implementation=generate(
	<<'END_OF_DERIVATIVE_FUNC_IMPLEMENTATION','','','');
function_options & function_options::derivative_func(derivative_funcp_${N} d)
{
	test_and_set_nparams(${N});
	derivative_f = derivative_funcp(d);
	return *this;
}
END_OF_DERIVATIVE_FUNC_IMPLEMENTATION

$power_func_implementation=generate(
	<<'END_OF_POWER_FUNC_IMPLEMENTATION','','','');
function_options & function_options::power_func(power_funcp_${N} d)
{
	test_and_set_nparams(${N});
	power_f = power_funcp(d);
	return *this;
}
END_OF_POWER_FUNC_IMPLEMENTATION

$series_func_implementation=generate(
	<<'END_OF_SERIES_FUNC_IMPLEMENTATION','','','');
function_options & function_options::series_func(series_funcp_${N} s)
{
	test_and_set_nparams(${N});
	series_f = series_funcp(s);
	return *this;
}
END_OF_SERIES_FUNC_IMPLEMENTATION

$interface=<<END_OF_INTERFACE;
/** \@file function.h
 *
 *  Interface to class of symbolic functions. */

/*
 *  This file was generated automatically by function.pl.
 *  Please do not modify it directly, edit the perl script instead!
 *  function.pl options: \$maxargs=${maxargs}
 *
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

#ifndef __GINAC_FUNCTION_H__
#define __GINAC_FUNCTION_H__

#include <string>
#include <vector>

// CINT needs <algorithm> to work properly with <vector>
#include <algorithm>

#include "exprseq.h"

// the following lines have been generated for max. ${maxargs} parameters
$declare_function_macro
// end of generated lines

#define REGISTER_FUNCTION(NAME,OPT) \\
unsigned NAME##_SERIAL::serial = \\
	GiNaC::function::register_new(GiNaC::function_options(#NAME, NAME##_NPARAMS).OPT);

namespace GiNaC {

class function;
class symmetry;

typedef ex (* eval_funcp)();
typedef ex (* evalf_funcp)();
typedef ex (* conjugate_funcp)();
typedef ex (* real_part_funcp)();
typedef ex (* imag_part_funcp)();
typedef ex (* derivative_funcp)();
typedef ex (* power_funcp)();
typedef ex (* series_funcp)();
typedef void (* print_funcp)();

// the following lines have been generated for max. ${maxargs} parameters
$typedef_eval_funcp
$typedef_evalf_funcp
$typedef_conjugate_funcp
$typedef_real_part_funcp
$typedef_imag_part_funcp
$typedef_derivative_funcp
$typedef_power_funcp
$typedef_series_funcp
$typedef_print_funcp
// end of generated lines

// Alternatively, an exvector may be passed into the static function, instead
// of individual ex objects.  Then, the number of arguments is not limited.
typedef ex (* eval_funcp_exvector)(const exvector &);
typedef ex (* evalf_funcp_exvector)(const exvector &);
typedef ex (* conjugate_funcp_exvector)(const exvector &);
typedef ex (* real_part_funcp_exvector)(const exvector &);
typedef ex (* imag_part_funcp_exvector)(const exvector &);
typedef ex (* derivative_funcp_exvector)(const exvector &, unsigned);
typedef ex (* power_funcp_exvector)(const exvector &, const ex &);
typedef ex (* series_funcp_exvector)(const exvector &, const relational &, int, unsigned);
typedef void (* print_funcp_exvector)(const exvector &, const print_context &);


class function_options
{
	friend class function;
	friend class fderivative;
public:
	function_options();
	function_options(std::string const & n, std::string const & tn=std::string());
	function_options(std::string const & n, unsigned np);
	~function_options();
	void initialize();

	function_options & dummy() { return *this; }
	function_options & set_name(std::string const & n, std::string const & tn=std::string());
	function_options & latex_name(std::string const & tn);
// the following lines have been generated for max. ${maxargs} parameters
$eval_func_interface
$evalf_func_interface
$conjugate_func_interface
$real_part_func_interface
$imag_part_func_interface
$derivative_func_interface
$power_func_interface
$series_func_interface
$print_func_interface
// end of generated lines
	function_options & eval_func(eval_funcp_exvector e);
	function_options & evalf_func(evalf_funcp_exvector ef);
	function_options & conjugate_func(conjugate_funcp_exvector d);
	function_options & real_part_func(real_part_funcp_exvector d);
	function_options & imag_part_func(imag_part_funcp_exvector d);
	function_options & derivative_func(derivative_funcp_exvector d);
	function_options & power_func(power_funcp_exvector d);
	function_options & series_func(series_funcp_exvector s);

	template <class Ctx> function_options & print_func(print_funcp_exvector p)
	{
		print_use_exvector_args = true;
		set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
		return *this;
	}

	function_options & set_return_type(unsigned rt, tinfo_t rtt=NULL);
	function_options & do_not_evalf_params();
	function_options & remember(unsigned size, unsigned assoc_size=0,
	                            unsigned strategy=remember_strategies::delete_never);
	function_options & overloaded(unsigned o);
	function_options & set_symmetry(const symmetry & s);

	std::string get_name() const { return name; }
	unsigned get_nparams() const { return nparams; }

protected:
	bool has_derivative() const { return derivative_f != NULL; }
	bool has_power() const { return power_f != NULL; }
	void test_and_set_nparams(unsigned n);
	void set_print_func(unsigned id, print_funcp f);

	std::string name;
	std::string TeX_name;

	unsigned nparams;

	eval_funcp eval_f;
	evalf_funcp evalf_f;
	conjugate_funcp conjugate_f;
	real_part_funcp real_part_f;
	imag_part_funcp imag_part_f;
	derivative_funcp derivative_f;
	power_funcp power_f;
	series_funcp series_f;
	std::vector<print_funcp> print_dispatch_table;

	bool evalf_params_first;

	bool use_return_type;
	unsigned return_type;
	tinfo_t return_type_tinfo;

	bool use_remember;
	unsigned remember_size;
	unsigned remember_assoc_size;
	unsigned remember_strategy;

	bool eval_use_exvector_args;
	bool evalf_use_exvector_args;
	bool conjugate_use_exvector_args;
	bool real_part_use_exvector_args;
	bool imag_part_use_exvector_args;
	bool derivative_use_exvector_args;
	bool power_use_exvector_args;
	bool series_use_exvector_args;
	bool print_use_exvector_args;

	unsigned functions_with_same_name;

	ex symtree;
};


/** Exception class thrown by classes which provide their own series expansion
 *  to signal that ordinary Taylor expansion is safe. */
class do_taylor {};


/** The class function is used to implement builtin functions like sin, cos...
	and user defined functions */
class function : public exprseq
{
	GINAC_DECLARE_REGISTERED_CLASS(function, exprseq)

	// CINT has a linking problem
#ifndef __MAKECINT__
	friend void ginsh_get_ginac_functions();
#endif // def __MAKECINT__

	friend class remember_table_entry;
	// friend class remember_table_list;
	// friend class remember_table;

// member functions

	// other constructors
public:
	function(unsigned ser);
	// the following lines have been generated for max. ${maxargs} parameters
$constructors_interface
	// end of generated lines
	function(unsigned ser, const exprseq & es);
	function(unsigned ser, const exvector & v, bool discardable = false);
	function(unsigned ser, std::auto_ptr<exvector> vp);

	// functions overriding virtual functions from base classes
public:
	void print(const print_context & c, unsigned level = 0) const;
	unsigned precedence() const {return 70;}
	ex expand(unsigned options=0) const;
	ex eval(int level=0) const;
	ex evalf(int level=0) const;
	unsigned calchash() const;
	ex series(const relational & r, int order, unsigned options = 0) const;
	ex thiscontainer(const exvector & v) const;
	ex thiscontainer(std::auto_ptr<exvector> vp) const;
	ex conjugate() const;
	ex real_part() const;
	ex imag_part() const;
protected:
	ex derivative(const symbol & s) const;
	bool is_equal_same_type(const basic & other) const;
	bool match_same_type(const basic & other) const;
	unsigned return_type() const;
	tinfo_t return_type_tinfo() const;
	
	// new virtual functions which can be overridden by derived classes
	// none
	
	// non-virtual functions in this class
protected:
	ex pderivative(unsigned diff_param) const; // partial differentiation
	static std::vector<function_options> & registered_functions();
	bool lookup_remember_table(ex & result) const;
	void store_remember_table(ex const & result) const;
public:
	ex power(const ex & exp) const;
	static unsigned register_new(function_options const & opt);
	static unsigned current_serial;
	static unsigned find_function(const std::string &name, unsigned nparams);
	unsigned get_serial() const {return serial;}
	std::string get_name() const;

// member variables

protected:
	unsigned serial;
};

// utility functions/macros

template <typename T>
inline bool is_the_function(const ex & x)
{
	return is_exactly_a<function>(x)
	    && ex_to<function>(x).get_serial() == T::serial;
}

// Check whether OBJ is the specified symbolic function.
#define is_ex_the_function(OBJ, FUNCNAME) (GiNaC::is_the_function<FUNCNAME##_SERIAL>(OBJ))

} // namespace GiNaC

#endif // ndef __GINAC_FUNCTION_H__

END_OF_INTERFACE

$implementation=<<END_OF_IMPLEMENTATION;
/** \@file function.cpp
 *
 *  Implementation of class of symbolic functions. */

/*
 *  This file was generated automatically by function.pl.
 *  Please do not modify it directly, edit the perl script instead!
 *  function.pl options: \$maxargs=${maxargs}
 *
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

#include <iostream>
#include <string>
#include <stdexcept>
#include <list>
#include <limits>

#include "function.h"
#include "operators.h"
#include "fderivative.h"
#include "ex.h"
#include "lst.h"
#include "symmetry.h"
#include "print.h"
#include "power.h"
#include "archive.h"
#include "inifcns.h"
#include "tostring.h"
#include "utils.h"
#include "remember.h"

namespace GiNaC {

//////////
// helper class function_options
//////////

function_options::function_options()
{
	initialize();
}

function_options::function_options(std::string const & n, std::string const & tn)
{
	initialize();
	set_name(n, tn);
}

function_options::function_options(std::string const & n, unsigned np)
{
	initialize();
	set_name(n, std::string());
	nparams = np;
}

function_options::~function_options()
{
	// nothing to clean up at the moment
}

void function_options::initialize()
{
	set_name("unnamed_function", "\\\\mbox{unnamed}");
	nparams = 0;
	eval_f = evalf_f = real_part_f = imag_part_f = conjugate_f = derivative_f
		= power_f = series_f = 0;
	evalf_params_first = true;
	use_return_type = false;
	eval_use_exvector_args = false;
	evalf_use_exvector_args = false;
	conjugate_use_exvector_args = false;
	real_part_use_exvector_args = false;
	imag_part_use_exvector_args = false;
	derivative_use_exvector_args = false;
	power_use_exvector_args = false;
	series_use_exvector_args = false;
	print_use_exvector_args = false;
	use_remember = false;
	functions_with_same_name = 1;
	symtree = 0;
}

function_options & function_options::set_name(std::string const & n,
                                              std::string const & tn)
{
	name = n;
	if (tn==std::string())
		TeX_name = "\\\\mbox{"+name+"}";
	else
		TeX_name = tn;
	return *this;
}

function_options & function_options::latex_name(std::string const & tn)
{
	TeX_name = tn;
	return *this;
}

// the following lines have been generated for max. ${maxargs} parameters
$eval_func_implementation
$evalf_func_implementation
$conjugate_func_implementation
$real_part_func_implementation
$imag_part_func_implementation
$derivative_func_implementation
$power_func_implementation
$series_func_implementation
// end of generated lines

function_options& function_options::eval_func(eval_funcp_exvector e)
{
	eval_use_exvector_args = true;
	eval_f = eval_funcp(e);
	return *this;
}
function_options& function_options::evalf_func(evalf_funcp_exvector ef)
{
	evalf_use_exvector_args = true;
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options& function_options::conjugate_func(conjugate_funcp_exvector c)
{
	conjugate_use_exvector_args = true;
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options& function_options::real_part_func(real_part_funcp_exvector c)
{
	real_part_use_exvector_args = true;
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options& function_options::imag_part_func(imag_part_funcp_exvector c)
{
	imag_part_use_exvector_args = true;
	imag_part_f = imag_part_funcp(c);
	return *this;
}

function_options& function_options::derivative_func(derivative_funcp_exvector d)
{
	derivative_use_exvector_args = true;
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options& function_options::power_func(power_funcp_exvector d)
{
	power_use_exvector_args = true;
	power_f = power_funcp(d);
	return *this;
}
function_options& function_options::series_func(series_funcp_exvector s)
{
	series_use_exvector_args = true;
	series_f = series_funcp(s);
	return *this;
}

function_options & function_options::set_return_type(unsigned rt, tinfo_t rtt)
{
	use_return_type = true;
	return_type = rt;
	return_type_tinfo = rtt;
	return *this;
}

function_options & function_options::do_not_evalf_params()
{
	evalf_params_first = false;
	return *this;
}

function_options & function_options::remember(unsigned size,
                                              unsigned assoc_size,
                                              unsigned strategy)
{
	use_remember = true;
	remember_size = size;
	remember_assoc_size = assoc_size;
	remember_strategy = strategy;
	return *this;
}

function_options & function_options::overloaded(unsigned o)
{
	functions_with_same_name = o;
	return *this;
}

function_options & function_options::set_symmetry(const symmetry & s)
{
	symtree = s;
	return *this;
}
	
void function_options::test_and_set_nparams(unsigned n)
{
	if (nparams==0) {
		nparams = n;
	} else if (nparams!=n) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// be caught anyhow
		std::cerr << "WARNING: " << name << "(): number of parameters ("
		          << n << ") differs from number set before (" 
		          << nparams << ")" << std::endl;
	}
}

void function_options::set_print_func(unsigned id, print_funcp f)
{
	if (id >= print_dispatch_table.size())
		print_dispatch_table.resize(id + 1);
	print_dispatch_table[id] = f;
}

/** This can be used as a hook for external applications. */
unsigned function::current_serial = 0;


GINAC_IMPLEMENT_REGISTERED_CLASS(function, exprseq)

//////////
// default constructor
//////////

// public

function::function() : serial(0)
{
	tinfo_key = &function::tinfo_static;
}

//////////
// other constructors
//////////

// public

function::function(unsigned ser) : serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

// the following lines have been generated for max. ${maxargs} parameters
$constructors_implementation
// end of generated lines

function::function(unsigned ser, const exprseq & es) : exprseq(es), serial(ser)
{
	tinfo_key = &function::tinfo_static;

	// Force re-evaluation even if the exprseq was already evaluated
	// (the exprseq copy constructor copies the flags)
	clearflag(status_flags::evaluated);
}

function::function(unsigned ser, const exvector & v, bool discardable) 
  : exprseq(v,discardable), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

function::function(unsigned ser, std::auto_ptr<exvector> vp) 
  : exprseq(vp), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

//////////
// archiving
//////////

/** Construct object from archive_node. */
function::function(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	// Find serial number by function name
	std::string s;
	if (n.find_string("name", s)) {
		unsigned int ser = 0;
		std::vector<function_options>::const_iterator i = registered_functions().begin(), iend = registered_functions().end();
		while (i != iend) {
			if (s == i->name) {
				serial = ser;
				return;
			}
			++i; ++ser;
		}
		throw (std::runtime_error("unknown function '" + s + "' in archive"));
	} else
		throw (std::runtime_error("unnamed function in archive"));
}

/** Unarchive the object. */
ex function::unarchive(const archive_node &n, lst &sym_lst)
{
	return (new function(n, sym_lst))->setflag(status_flags::dynallocated);
}

/** Archive the object. */
void function::archive(archive_node &n) const
{
	inherited::archive(n);
	GINAC_ASSERT(serial < registered_functions().size());
	n.add_string("name", registered_functions()[serial].name);
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void function::print(const print_context & c, unsigned level) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	const std::vector<print_funcp> &pdt = opt.print_dispatch_table;

	// Dynamically dispatch on print_context type
	const print_context_class_info *pc_info = &c.get_class_info();

next_context:
	unsigned id = pc_info->options.get_id();
	if (id >= pdt.size() || pdt[id] == NULL) {

		// Method not found, try parent print_context class
		const print_context_class_info *parent_pc_info = pc_info->get_parent();
		if (parent_pc_info) {
			pc_info = parent_pc_info;
			goto next_context;
		}

		// Method still not found, use default output
		if (is_a<print_tree>(c)) {

			c.s << std::string(level, ' ') << class_name() << " "
			    << opt.name << " @" << this
			    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
			    << ", nops=" << nops()
			    << std::endl;
			unsigned delta_indent = static_cast<const print_tree &>(c).delta_indent;
			for (size_t i=0; i<seq.size(); ++i)
				seq[i].print(c, level + delta_indent);
			c.s << std::string(level + delta_indent, ' ') << "=====" << std::endl;

		} else if (is_a<print_csrc>(c)) {

			// Print function name in lowercase
			std::string lname = opt.name;
			size_t num = lname.size();
			for (size_t i=0; i<num; i++)
				lname[i] = tolower(lname[i]);
			c.s << lname;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());

		} else if (is_a<print_latex>(c)) {
			c.s << opt.TeX_name;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
		} else {
			c.s << opt.name;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
		}

	} else {

		// Method found, call it
		current_serial = serial;
		if (opt.print_use_exvector_args)
			((print_funcp_exvector)pdt[id])(seq, c);
		else switch (opt.nparams) {
			// the following lines have been generated for max. ${maxargs} parameters
${print_switch_statement}
			// end of generated lines
		default:
			throw(std::logic_error("function::print(): invalid nparams"));
		}
	}
}

ex function::expand(unsigned options) const
{
	// Only expand arguments when asked to do so
	if (options & expand_options::expand_function_args)
		return inherited::expand(options);
	else
		return (options == 0) ? setflag(status_flags::expanded) : *this;
}

ex function::eval(int level) const
{
	if (level>1) {
		// first evaluate children, then we will end up here again
		return function(serial,evalchildren(level));
	}

	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	// Canonicalize argument order according to the symmetry properties
	if (seq.size() > 1 && !(opt.symtree.is_zero())) {
		exvector v = seq;
		GINAC_ASSERT(is_a<symmetry>(opt.symtree));
		int sig = canonicalize(v.begin(), ex_to<symmetry>(opt.symtree));
		if (sig != std::numeric_limits<int>::max()) {
			// Something has changed while sorting arguments, more evaluations later
			if (sig == 0)
				return _ex0;
			return ex(sig) * thiscontainer(v);
		}
	}

	if (opt.eval_f==0) {
		return this->hold();
	}

	bool use_remember = opt.use_remember;
	ex eval_result;
	if (use_remember && lookup_remember_table(eval_result)) {
		return eval_result;
	}
	current_serial = serial;
	if (opt.eval_use_exvector_args)
		eval_result = ((eval_funcp_exvector)(opt.eval_f))(seq);
	else
	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${eval_switch_statement}
		// end of generated lines
	default:
		throw(std::logic_error("function::eval(): invalid nparams"));
	}
	if (use_remember) {
		store_remember_table(eval_result);
	}
	return eval_result;
}

ex function::evalf(int level) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	// Evaluate children first
	exvector eseq;
	if (level == 1 || !(opt.evalf_params_first))
		eseq = seq;
	else if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	else {
		eseq.reserve(seq.size());
		--level;
		exvector::const_iterator it = seq.begin(), itend = seq.end();
		while (it != itend) {
			eseq.push_back(it->evalf(level));
			++it;
		}
	}

	if (opt.evalf_f==0) {
		return function(serial,eseq).hold();
	}
	current_serial = serial;
	if (opt.evalf_use_exvector_args)
		return ((evalf_funcp_exvector)(opt.evalf_f))(seq);
	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${evalf_switch_statement}
		// end of generated lines
	}
	throw(std::logic_error("function::evalf(): invalid nparams"));
}

unsigned function::calchash() const
{
	unsigned v = golden_ratio_hash(golden_ratio_hash((p_int)tinfo()) ^ serial);
	for (size_t i=0; i<nops(); i++) {
		v = rotate_left(v);
		v ^= this->op(i).gethash();
	}

	if (flags & status_flags::evaluated) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}
	return v;
}

ex function::thiscontainer(const exvector & v) const
{
	return function(serial, v);
}

ex function::thiscontainer(std::auto_ptr<exvector> vp) const
{
	return function(serial, vp);
}

/** Implementation of ex::series for functions.
 *  \@see ex::series */
ex function::series(const relational & r, int order, unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.series_f==0) {
		return basic::series(r, order);
	}
	ex res;
	current_serial = serial;
	if (opt.series_use_exvector_args) {
		try {
			res = ((series_funcp_exvector)(opt.series_f))(seq, r, order, options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	}
	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${series_switch_statement}
		// end of generated lines
	}
	throw(std::logic_error("function::series(): invalid nparams"));
}

/** Implementation of ex::conjugate for functions. */
ex function::conjugate() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.conjugate_f==0) {
		return exprseq::conjugate();
	}

	if (opt.conjugate_use_exvector_args) {
		return ((conjugate_funcp_exvector)(opt.conjugate_f))(seq);
	}

	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${conjugate_switch_statement}
		// end of generated lines
	}
	throw(std::logic_error("function::conjugate(): invalid nparams"));
}

/** Implementation of ex::real_part for functions. */
ex function::real_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.real_part_f==0)
		return basic::real_part();

	if (opt.real_part_use_exvector_args)
		return ((real_part_funcp_exvector)(opt.real_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${real_part_switch_statement}
		// end of generated lines
	}
	throw(std::logic_error("function::real_part(): invalid nparams"));
}

/** Implementation of ex::imag_part for functions. */
ex function::imag_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.imag_part_f==0)
		return basic::imag_part();

	if (opt.imag_part_use_exvector_args)
		return ((imag_part_funcp_exvector)(opt.imag_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${imag_part_switch_statement}
		// end of generated lines
	}
	throw(std::logic_error("function::imag_part(): invalid nparams"));
}

// protected

/** Implementation of ex::diff() for functions. It applies the chain rule,
 *  except for the Order term function.
 *  \@see ex::diff */
ex function::derivative(const symbol & s) const
{
	ex result;

	if (serial == Order_SERIAL::serial) {
		// Order Term function only differentiates the argument
		return Order(seq[0].diff(s));
	} else {
		// Chain rule
		ex arg_diff;
		size_t num = seq.size();
		for (size_t i=0; i<num; i++) {
			arg_diff = seq[i].diff(s);
			// We apply the chain rule only when it makes sense.  This is not
			// just for performance reasons but also to allow functions to
			// throw when differentiated with respect to one of its arguments
			// without running into trouble with our automatic full
			// differentiation:
			if (!arg_diff.is_zero())
				result += pderivative(i)*arg_diff;
		}
	}
	return result;
}

int function::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return serial < o.serial ? -1 : 1;
	else
		return exprseq::compare_same_type(o);
}

bool function::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return false;
	else
		return exprseq::is_equal_same_type(o);
}

bool function::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	return serial == o.serial;
}

unsigned function::return_type() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type;
	} else {
		// Default behavior is to use the return type of the first
		// argument. Thus, exp() of a matrix behaves like a matrix, etc.
		if (seq.empty())
			return return_types::commutative;
		else
			return seq.begin()->return_type();
	}
}

tinfo_t function::return_type_tinfo() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type_tinfo;
	} else {
		// Default behavior is to use the return type of the first
		// argument. Thus, exp() of a matrix behaves like a matrix, etc.
		if (seq.empty())
			return this;
		else
			return seq.begin()->return_type_tinfo();
	}
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// protected

ex function::pderivative(unsigned diff_param) const // partial differentiation
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	// No derivative defined? Then return abstract derivative object
	if (opt.derivative_f == NULL)
		return fderivative(serial, diff_param, seq);

	current_serial = serial;
	if (opt.derivative_use_exvector_args)
		return ((derivative_funcp_exvector)(opt.derivative_f))(seq, diff_param);
	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${diff_switch_statement}
		// end of generated lines
	}
	throw(std::logic_error("function::pderivative(): no diff function defined"));
}

ex function::power(const ex & power_param) const // power of function
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	// No derivative defined? Then return abstract derivative object
	if (opt.power_f == NULL)
		return (new power::power(*this, power_param))->setflag(status_flags::dynallocated |
	                                               status_flags::evaluated);

	current_serial = serial;
	if (opt.power_use_exvector_args)
		return ((power_funcp_exvector)(opt.power_f))(seq,  power_param);
	switch (opt.nparams) {
		// the following lines have been generated for max. ${maxargs} parameters
${power_switch_statement}
		// end of generated lines
	}
	throw(std::logic_error("function::power(): no power function defined"));
}

std::vector<function_options> & function::registered_functions()
{
	static std::vector<function_options> * rf = new std::vector<function_options>;
	return *rf;
}

bool function::lookup_remember_table(ex & result) const
{
	return remember_table::remember_tables()[this->serial].lookup_entry(*this,result);
}

void function::store_remember_table(ex const & result) const
{
	remember_table::remember_tables()[this->serial].add_entry(*this,result);
}

// public

unsigned function::register_new(function_options const & opt)
{
	size_t same_name = 0;
	for (size_t i=0; i<registered_functions().size(); ++i) {
		if (registered_functions()[i].name==opt.name) {
			++same_name;
		}
	}
	if (same_name>=opt.functions_with_same_name) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// caught anyhow
		std::cerr << "WARNING: function name " << opt.name
		          << " already in use!" << std::endl;
	}
	registered_functions().push_back(opt);
	if (opt.use_remember) {
		remember_table::remember_tables().
			push_back(remember_table(opt.remember_size,
			                         opt.remember_assoc_size,
			                         opt.remember_strategy));
	} else {
		remember_table::remember_tables().push_back(remember_table());
	}
	return registered_functions().size()-1;
}

/** Find serial number of function by name and number of parameters.
 *  Throws exception if function was not found. */
unsigned function::find_function(const std::string &name, unsigned nparams)
{
	std::vector<function_options>::const_iterator i = function::registered_functions().begin(), end = function::registered_functions().end();
	unsigned serial = 0;
	while (i != end) {
		if (i->get_name() == name && i->get_nparams() == nparams)
			return serial;
		++i;
		++serial;
	}
	throw (std::runtime_error("no function '" + name + "' with " + ToString(nparams) + " parameters defined"));
}

/** Return the print name of the function. */
std::string function::get_name() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	return registered_functions()[serial].name;
}

} // namespace GiNaC

END_OF_IMPLEMENTATION

print "Creating interface file function.h...";
open OUT,">function.h" or die "cannot open function.h";
print OUT $interface;
close OUT;
print "ok.\n";

print "Creating implementation file function.cpp...";
open OUT,">function.cpp" or die "cannot open function.cpp";
print OUT $implementation;
close OUT;
print "ok.\n";

print "done.\n";
