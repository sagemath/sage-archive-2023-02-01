/** @file function.cpp
 *
 *  Implementation of class of symbolic functions. */

/*
 *  This file was generated automatically by function.pl.
 *  Please do not modify it directly, edit the perl script instead!
 *  function.pl options: $maxargs=14
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
#include "relational.h"
#include "archive.h"
#include "inifcns.h"
#include "tostring.h"
#include "utils.h"
#include "remember.h"

extern "C" {
	PyObject* exvector_to_PyTuple(GiNaC::exvector seq);
	GiNaC::ex pyExpression_to_ex(PyObject* s);
	PyObject* ex_to_pyExpression(GiNaC::ex e);
}
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
	set_name("unnamed_function", "\\mbox{unnamed}");
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
	python_func = false;
	functions_with_same_name = 1;
	symtree = 0;
}

function_options & function_options::set_name(std::string const & n,
                                              std::string const & tn)
{
	name = n;
	if (tn==std::string())
		TeX_name = "\\mbox{"+name+"}";
	else
		TeX_name = tn;
	return *this;
}

function_options & function_options::latex_name(std::string const & tn)
{
	TeX_name = tn;
	return *this;
}

// the following lines have been generated for max. 14 parameters
function_options & function_options::eval_func(eval_funcp_1 e)
{
	test_and_set_nparams(1);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_2 e)
{
	test_and_set_nparams(2);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_3 e)
{
	test_and_set_nparams(3);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_4 e)
{
	test_and_set_nparams(4);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_5 e)
{
	test_and_set_nparams(5);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_6 e)
{
	test_and_set_nparams(6);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_7 e)
{
	test_and_set_nparams(7);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_8 e)
{
	test_and_set_nparams(8);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_9 e)
{
	test_and_set_nparams(9);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_10 e)
{
	test_and_set_nparams(10);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_11 e)
{
	test_and_set_nparams(11);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_12 e)
{
	test_and_set_nparams(12);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_13 e)
{
	test_and_set_nparams(13);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_14 e)
{
	test_and_set_nparams(14);
	eval_f = eval_funcp(e);
	return *this;
}

function_options & function_options::evalf_func(evalf_funcp_1 ef)
{
	test_and_set_nparams(1);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_2 ef)
{
	test_and_set_nparams(2);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_3 ef)
{
	test_and_set_nparams(3);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_4 ef)
{
	test_and_set_nparams(4);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_5 ef)
{
	test_and_set_nparams(5);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_6 ef)
{
	test_and_set_nparams(6);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_7 ef)
{
	test_and_set_nparams(7);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_8 ef)
{
	test_and_set_nparams(8);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_9 ef)
{
	test_and_set_nparams(9);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_10 ef)
{
	test_and_set_nparams(10);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_11 ef)
{
	test_and_set_nparams(11);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_12 ef)
{
	test_and_set_nparams(12);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_13 ef)
{
	test_and_set_nparams(13);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_14 ef)
{
	test_and_set_nparams(14);
	evalf_f = evalf_funcp(ef);
	return *this;
}

function_options & function_options::conjugate_func(conjugate_funcp_1 c)
{
	test_and_set_nparams(1);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_2 c)
{
	test_and_set_nparams(2);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_3 c)
{
	test_and_set_nparams(3);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_4 c)
{
	test_and_set_nparams(4);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_5 c)
{
	test_and_set_nparams(5);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_6 c)
{
	test_and_set_nparams(6);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_7 c)
{
	test_and_set_nparams(7);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_8 c)
{
	test_and_set_nparams(8);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_9 c)
{
	test_and_set_nparams(9);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_10 c)
{
	test_and_set_nparams(10);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_11 c)
{
	test_and_set_nparams(11);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_12 c)
{
	test_and_set_nparams(12);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_13 c)
{
	test_and_set_nparams(13);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_14 c)
{
	test_and_set_nparams(14);
	conjugate_f = conjugate_funcp(c);
	return *this;
}

function_options & function_options::real_part_func(real_part_funcp_1 c)
{
	test_and_set_nparams(1);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_2 c)
{
	test_and_set_nparams(2);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_3 c)
{
	test_and_set_nparams(3);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_4 c)
{
	test_and_set_nparams(4);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_5 c)
{
	test_and_set_nparams(5);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_6 c)
{
	test_and_set_nparams(6);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_7 c)
{
	test_and_set_nparams(7);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_8 c)
{
	test_and_set_nparams(8);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_9 c)
{
	test_and_set_nparams(9);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_10 c)
{
	test_and_set_nparams(10);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_11 c)
{
	test_and_set_nparams(11);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_12 c)
{
	test_and_set_nparams(12);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_13 c)
{
	test_and_set_nparams(13);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_14 c)
{
	test_and_set_nparams(14);
	real_part_f = real_part_funcp(c);
	return *this;
}

function_options & function_options::imag_part_func(imag_part_funcp_1 c)
{
	test_and_set_nparams(1);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_2 c)
{
	test_and_set_nparams(2);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_3 c)
{
	test_and_set_nparams(3);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_4 c)
{
	test_and_set_nparams(4);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_5 c)
{
	test_and_set_nparams(5);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_6 c)
{
	test_and_set_nparams(6);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_7 c)
{
	test_and_set_nparams(7);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_8 c)
{
	test_and_set_nparams(8);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_9 c)
{
	test_and_set_nparams(9);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_10 c)
{
	test_and_set_nparams(10);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_11 c)
{
	test_and_set_nparams(11);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_12 c)
{
	test_and_set_nparams(12);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_13 c)
{
	test_and_set_nparams(13);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_14 c)
{
	test_and_set_nparams(14);
	imag_part_f = imag_part_funcp(c);
	return *this;
}

function_options & function_options::derivative_func(derivative_funcp_1 d)
{
	test_and_set_nparams(1);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_2 d)
{
	test_and_set_nparams(2);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_3 d)
{
	test_and_set_nparams(3);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_4 d)
{
	test_and_set_nparams(4);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_5 d)
{
	test_and_set_nparams(5);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_6 d)
{
	test_and_set_nparams(6);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_7 d)
{
	test_and_set_nparams(7);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_8 d)
{
	test_and_set_nparams(8);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_9 d)
{
	test_and_set_nparams(9);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_10 d)
{
	test_and_set_nparams(10);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_11 d)
{
	test_and_set_nparams(11);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_12 d)
{
	test_and_set_nparams(12);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_13 d)
{
	test_and_set_nparams(13);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_14 d)
{
	test_and_set_nparams(14);
	derivative_f = derivative_funcp(d);
	return *this;
}

function_options & function_options::power_func(power_funcp_1 d)
{
	test_and_set_nparams(1);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_2 d)
{
	test_and_set_nparams(2);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_3 d)
{
	test_and_set_nparams(3);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_4 d)
{
	test_and_set_nparams(4);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_5 d)
{
	test_and_set_nparams(5);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_6 d)
{
	test_and_set_nparams(6);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_7 d)
{
	test_and_set_nparams(7);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_8 d)
{
	test_and_set_nparams(8);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_9 d)
{
	test_and_set_nparams(9);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_10 d)
{
	test_and_set_nparams(10);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_11 d)
{
	test_and_set_nparams(11);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_12 d)
{
	test_and_set_nparams(12);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_13 d)
{
	test_and_set_nparams(13);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_14 d)
{
	test_and_set_nparams(14);
	power_f = power_funcp(d);
	return *this;
}

function_options & function_options::series_func(series_funcp_1 s)
{
	test_and_set_nparams(1);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_2 s)
{
	test_and_set_nparams(2);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_3 s)
{
	test_and_set_nparams(3);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_4 s)
{
	test_and_set_nparams(4);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_5 s)
{
	test_and_set_nparams(5);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_6 s)
{
	test_and_set_nparams(6);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_7 s)
{
	test_and_set_nparams(7);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_8 s)
{
	test_and_set_nparams(8);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_9 s)
{
	test_and_set_nparams(9);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_10 s)
{
	test_and_set_nparams(10);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_11 s)
{
	test_and_set_nparams(11);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_12 s)
{
	test_and_set_nparams(12);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_13 s)
{
	test_and_set_nparams(13);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_14 s)
{
	test_and_set_nparams(14);
	series_f = series_funcp(s);
	return *this;
}

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

function_options& function_options::eval_func(PyObject* e)
{
	eval_f = eval_funcp(e);
	return *this;
}
function_options& function_options::evalf_func(PyObject* ef)
{
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options& function_options::conjugate_func(PyObject* c)
{
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options& function_options::real_part_func(PyObject* c)
{
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options& function_options::imag_part_func(PyObject* c)
{
	imag_part_f = imag_part_funcp(c);
	return *this;
}

function_options& function_options::derivative_func(PyObject* d)
{
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options& function_options::power_func(PyObject* d)
{
	power_f = power_funcp(d);
	return *this;
}
function_options& function_options::series_func(PyObject* s)
{
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

// the following lines have been generated for max. 14 parameters
function::function(unsigned ser, const ex & param1)
	: exprseq(param1), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2)
	: exprseq(param1, param2), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3)
	: exprseq(param1, param2, param3), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4)
	: exprseq(param1, param2, param3, param4), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5)
	: exprseq(param1, param2, param3, param4, param5), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6)
	: exprseq(param1, param2, param3, param4, param5, param6), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7)
	: exprseq(param1, param2, param3, param4, param5, param6, param7), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8)
	: exprseq(param1, param2, param3, param4, param5, param6, param7, param8), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9)
	: exprseq(param1, param2, param3, param4, param5, param6, param7, param8, param9), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10)
	: exprseq(param1, param2, param3, param4, param5, param6, param7, param8, param9, param10), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11)
	: exprseq(param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12)
	: exprseq(param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12, const ex & param13)
	: exprseq(param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12, param13), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12, const ex & param13, const ex & param14)
	: exprseq(param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12, param13, param14), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

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
			// the following lines have been generated for max. 14 parameters
		case 1:
			((print_funcp_1)(pdt[id]))(seq[1-1], c);
			break;
		case 2:
			((print_funcp_2)(pdt[id]))(seq[1-1], seq[2-1], c);
			break;
		case 3:
			((print_funcp_3)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], c);
			break;
		case 4:
			((print_funcp_4)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], c);
			break;
		case 5:
			((print_funcp_5)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], c);
			break;
		case 6:
			((print_funcp_6)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], c);
			break;
		case 7:
			((print_funcp_7)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], c);
			break;
		case 8:
			((print_funcp_8)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], c);
			break;
		case 9:
			((print_funcp_9)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], c);
			break;
		case 10:
			((print_funcp_10)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], c);
			break;
		case 11:
			((print_funcp_11)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], c);
			break;
		case 12:
			((print_funcp_12)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], c);
			break;
		case 13:
			((print_funcp_13)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], c);
			break;
		case 14:
			((print_funcp_14)(pdt[id]))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1], c);
			break;

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

	if (opt.python_func && PyCallable_Check((PyObject*)opt.eval_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// call opt.eval_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.eval_f, 
				args, NULL);
		Py_DECREF(args);
		if (!pyresult) { 
			throw(std::runtime_error("function::eval(): python function raised exception"));
		}
		// convert output Expression to an ex
		eval_result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::eval(): python function (Expression_to_ex) raised exception"));
		}
	}
	else if (opt.eval_use_exvector_args)
		eval_result = ((eval_funcp_exvector)(opt.eval_f))(seq);
	else
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		eval_result = ((eval_funcp_1)(opt.eval_f))(seq[1-1]);
		break;
	case 2:
		eval_result = ((eval_funcp_2)(opt.eval_f))(seq[1-1], seq[2-1]);
		break;
	case 3:
		eval_result = ((eval_funcp_3)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1]);
		break;
	case 4:
		eval_result = ((eval_funcp_4)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1]);
		break;
	case 5:
		eval_result = ((eval_funcp_5)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1]);
		break;
	case 6:
		eval_result = ((eval_funcp_6)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1]);
		break;
	case 7:
		eval_result = ((eval_funcp_7)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1]);
		break;
	case 8:
		eval_result = ((eval_funcp_8)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1]);
		break;
	case 9:
		eval_result = ((eval_funcp_9)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1]);
		break;
	case 10:
		eval_result = ((eval_funcp_10)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1]);
		break;
	case 11:
		eval_result = ((eval_funcp_11)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1]);
		break;
	case 12:
		eval_result = ((eval_funcp_12)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1]);
		break;
	case 13:
		eval_result = ((eval_funcp_13)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1]);
		break;
	case 14:
		eval_result = ((eval_funcp_14)(opt.eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1]);
		break;

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
	if (opt.python_func && PyCallable_Check((PyObject*)opt.evalf_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// call opt.evalf_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.evalf_f, 
				args, NULL);
		Py_DECREF(args);
		if (!pyresult) { 
			throw(std::runtime_error("function::evalf(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::evalf(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.evalf_use_exvector_args)
		return ((evalf_funcp_exvector)(opt.evalf_f))(seq);
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return ((evalf_funcp_1)(opt.evalf_f))(eseq[1-1]);
	case 2:
		return ((evalf_funcp_2)(opt.evalf_f))(eseq[1-1], eseq[2-1]);
	case 3:
		return ((evalf_funcp_3)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1]);
	case 4:
		return ((evalf_funcp_4)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1]);
	case 5:
		return ((evalf_funcp_5)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1]);
	case 6:
		return ((evalf_funcp_6)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1]);
	case 7:
		return ((evalf_funcp_7)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1]);
	case 8:
		return ((evalf_funcp_8)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1], eseq[8-1]);
	case 9:
		return ((evalf_funcp_9)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1], eseq[8-1], eseq[9-1]);
	case 10:
		return ((evalf_funcp_10)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1], eseq[8-1], eseq[9-1], eseq[10-1]);
	case 11:
		return ((evalf_funcp_11)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1], eseq[8-1], eseq[9-1], eseq[10-1], eseq[11-1]);
	case 12:
		return ((evalf_funcp_12)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1], eseq[8-1], eseq[9-1], eseq[10-1], eseq[11-1], eseq[12-1]);
	case 13:
		return ((evalf_funcp_13)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1], eseq[8-1], eseq[9-1], eseq[10-1], eseq[11-1], eseq[12-1], eseq[13-1]);
	case 14:
		return ((evalf_funcp_14)(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], eseq[7-1], eseq[8-1], eseq[9-1], eseq[10-1], eseq[11-1], eseq[12-1], eseq[13-1], eseq[14-1]);

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
 *  @see ex::series */
ex function::series(const relational & r, int order, unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.series_f==0) {
		return basic::series(r, order);
	}
	ex res;
	current_serial = serial;
	if (opt.python_func && PyCallable_Check((PyObject*)opt.series_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// create a dictionary {'order': order, 'options':options}
		PyObject* kwds = Py_BuildValue("{s:i,s:I}","order",order,"options",options);
		// add variable to expand for as a keyword argument
		PyDict_SetItemString(kwds, "var", ex_to_pyExpression(r.lhs()));
		// add the point of expansion as a keyword argument
		PyDict_SetItemString(kwds, "at", ex_to_pyExpression(r.rhs()));
		// call opt.series_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.series_f, 
				args, kwds);
		Py_DECREF(args);
		Py_DECREF(kwds);
		if (!pyresult) { 
			throw(std::runtime_error("function::series(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::series(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.series_use_exvector_args) {
		try {
			res = ((series_funcp_exvector)(opt.series_f))(seq, r, order, options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	}
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		try {
			res = ((series_funcp_1)(opt.series_f))(seq[1-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 2:
		try {
			res = ((series_funcp_2)(opt.series_f))(seq[1-1], seq[2-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 3:
		try {
			res = ((series_funcp_3)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 4:
		try {
			res = ((series_funcp_4)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 5:
		try {
			res = ((series_funcp_5)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 6:
		try {
			res = ((series_funcp_6)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 7:
		try {
			res = ((series_funcp_7)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 8:
		try {
			res = ((series_funcp_8)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 9:
		try {
			res = ((series_funcp_9)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 10:
		try {
			res = ((series_funcp_10)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 11:
		try {
			res = ((series_funcp_11)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 12:
		try {
			res = ((series_funcp_12)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 13:
		try {
			res = ((series_funcp_13)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 14:
		try {
			res = ((series_funcp_14)(opt.series_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;

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

	if (opt.python_func && PyCallable_Check((PyObject*)opt.conjugate_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// call opt.conjugate_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.conjugate_f, 
				args, NULL);
		Py_DECREF(args);
		if (!pyresult) { 
			throw(std::runtime_error("function::conjugate(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::conjugate(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.conjugate_use_exvector_args) {
		return ((conjugate_funcp_exvector)(opt.conjugate_f))(seq);
	}

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return ((conjugate_funcp_1)(opt.conjugate_f))(seq[1-1]);
	case 2:
		return ((conjugate_funcp_2)(opt.conjugate_f))(seq[1-1], seq[2-1]);
	case 3:
		return ((conjugate_funcp_3)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1]);
	case 4:
		return ((conjugate_funcp_4)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1]);
	case 5:
		return ((conjugate_funcp_5)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1]);
	case 6:
		return ((conjugate_funcp_6)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1]);
	case 7:
		return ((conjugate_funcp_7)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1]);
	case 8:
		return ((conjugate_funcp_8)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1]);
	case 9:
		return ((conjugate_funcp_9)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1]);
	case 10:
		return ((conjugate_funcp_10)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1]);
	case 11:
		return ((conjugate_funcp_11)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1]);
	case 12:
		return ((conjugate_funcp_12)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1]);
	case 13:
		return ((conjugate_funcp_13)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1]);
	case 14:
		return ((conjugate_funcp_14)(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1]);

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

	if (opt.python_func && PyCallable_Check((PyObject*)opt.real_part_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// call opt.real_part_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.real_part_f, 
				args, NULL);
		Py_DECREF(args);
		if (!pyresult) { 
			throw(std::runtime_error("function::real_part(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::real_part(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.real_part_use_exvector_args)
		return ((real_part_funcp_exvector)(opt.real_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return ((real_part_funcp_1)(opt.real_part_f))(seq[1-1]);
	case 2:
		return ((real_part_funcp_2)(opt.real_part_f))(seq[1-1], seq[2-1]);
	case 3:
		return ((real_part_funcp_3)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1]);
	case 4:
		return ((real_part_funcp_4)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1]);
	case 5:
		return ((real_part_funcp_5)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1]);
	case 6:
		return ((real_part_funcp_6)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1]);
	case 7:
		return ((real_part_funcp_7)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1]);
	case 8:
		return ((real_part_funcp_8)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1]);
	case 9:
		return ((real_part_funcp_9)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1]);
	case 10:
		return ((real_part_funcp_10)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1]);
	case 11:
		return ((real_part_funcp_11)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1]);
	case 12:
		return ((real_part_funcp_12)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1]);
	case 13:
		return ((real_part_funcp_13)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1]);
	case 14:
		return ((real_part_funcp_14)(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1]);

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

	if (opt.python_func && PyCallable_Check((PyObject*)opt.imag_part_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// call opt.imag_part_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.imag_part_f, 
				args, NULL);
		Py_DECREF(args);
		if (!pyresult) { 
			throw(std::runtime_error("function::imag_part(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::imag_part(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.imag_part_use_exvector_args)
		return ((imag_part_funcp_exvector)(opt.imag_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return ((imag_part_funcp_1)(opt.imag_part_f))(seq[1-1]);
	case 2:
		return ((imag_part_funcp_2)(opt.imag_part_f))(seq[1-1], seq[2-1]);
	case 3:
		return ((imag_part_funcp_3)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1]);
	case 4:
		return ((imag_part_funcp_4)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1]);
	case 5:
		return ((imag_part_funcp_5)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1]);
	case 6:
		return ((imag_part_funcp_6)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1]);
	case 7:
		return ((imag_part_funcp_7)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1]);
	case 8:
		return ((imag_part_funcp_8)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1]);
	case 9:
		return ((imag_part_funcp_9)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1]);
	case 10:
		return ((imag_part_funcp_10)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1]);
	case 11:
		return ((imag_part_funcp_11)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1]);
	case 12:
		return ((imag_part_funcp_12)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1]);
	case 13:
		return ((imag_part_funcp_13)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1]);
	case 14:
		return ((imag_part_funcp_14)(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1]);

		// end of generated lines
	}
	throw(std::logic_error("function::imag_part(): invalid nparams"));
}

// protected

/** Implementation of ex::diff() for functions. It applies the chain rule,
 *  except for the Order term function.
 *  @see ex::diff */
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
	if (opt.python_func && PyCallable_Check((PyObject*)opt.derivative_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// create a dictionary {'diff_param': diff_param}
		PyObject* kwds = Py_BuildValue("{s:I}","diff_param",diff_param);
		// call opt.derivative_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.derivative_f, 
				args, kwds);
		Py_DECREF(args);
		Py_DECREF(kwds);
		if (!pyresult) { 
			throw(std::runtime_error("function::pderivative(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::pderivative(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.derivative_use_exvector_args)
		return ((derivative_funcp_exvector)(opt.derivative_f))(seq, diff_param);
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return ((derivative_funcp_1)(opt.derivative_f))(seq[1-1],diff_param);
	case 2:
		return ((derivative_funcp_2)(opt.derivative_f))(seq[1-1], seq[2-1],diff_param);
	case 3:
		return ((derivative_funcp_3)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1],diff_param);
	case 4:
		return ((derivative_funcp_4)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1],diff_param);
	case 5:
		return ((derivative_funcp_5)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1],diff_param);
	case 6:
		return ((derivative_funcp_6)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1],diff_param);
	case 7:
		return ((derivative_funcp_7)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1],diff_param);
	case 8:
		return ((derivative_funcp_8)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1],diff_param);
	case 9:
		return ((derivative_funcp_9)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1],diff_param);
	case 10:
		return ((derivative_funcp_10)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1],diff_param);
	case 11:
		return ((derivative_funcp_11)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1],diff_param);
	case 12:
		return ((derivative_funcp_12)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1],diff_param);
	case 13:
		return ((derivative_funcp_13)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1],diff_param);
	case 14:
		return ((derivative_funcp_14)(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1],diff_param);

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
	if (opt.python_func && PyCallable_Check((PyObject*)opt.power_f)) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = exvector_to_PyTuple(seq);
		// create a dictionary {'power_param': power_param}
		PyObject* kwds = PyDict_New();
		PyDict_SetItemString(kwds, "power_param", ex_to_pyExpression(power_param));
		// call opt.power_f with this list
		PyObject* pyresult = PyObject_Call((PyObject*)opt.power_f, 
				args, kwds);
		Py_DECREF(args);
		Py_DECREF(kwds);
		if (!pyresult) { 
			throw(std::runtime_error("function::power(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred()) { 
			throw(std::runtime_error("function::power(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.power_use_exvector_args)
		return ((power_funcp_exvector)(opt.power_f))(seq,  power_param);
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return ((power_funcp_1)(opt.power_f))(seq[1-1],power_param);
	case 2:
		return ((power_funcp_2)(opt.power_f))(seq[1-1], seq[2-1],power_param);
	case 3:
		return ((power_funcp_3)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1],power_param);
	case 4:
		return ((power_funcp_4)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1],power_param);
	case 5:
		return ((power_funcp_5)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1],power_param);
	case 6:
		return ((power_funcp_6)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1],power_param);
	case 7:
		return ((power_funcp_7)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1],power_param);
	case 8:
		return ((power_funcp_8)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1],power_param);
	case 9:
		return ((power_funcp_9)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1],power_param);
	case 10:
		return ((power_funcp_10)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1],power_param);
	case 11:
		return ((power_funcp_11)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1],power_param);
	case 12:
		return ((power_funcp_12)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1],power_param);
	case 13:
		return ((power_funcp_13)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1],power_param);
	case 14:
		return ((power_funcp_14)(opt.power_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], seq[7-1], seq[8-1], seq[9-1], seq[10-1], seq[11-1], seq[12-1], seq[13-1], seq[14-1],power_param);

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
		//
		// SAGE note: 
		// We suppress this warning since we allow a user to create 
		// functions with same name, but different number of args
		// Sage SFunction class checks existence of a function before
		// allocating a new one.
		//std::cerr << "WARNING: function name " << opt.name
		//          << " already in use!" << std::endl;
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

