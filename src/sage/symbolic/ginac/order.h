/** @file order.h
 *
 *  Definitions of order used for printing. */

/*
 *   Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
 *   Copyright (C) 2011 Jean-Pierre Flori <flori@enst.fr>
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

#ifndef __GINAC_ORDER_H__
#define __GINAC_ORDER_H__

#include "ex.h"
#include "basic.h"
#include "mul.h"
#include "power.h"
#include "add.h"
#include "symbol.h"
#include "function.h"
#include "fderivative.h"

namespace GiNaC {

struct ex_is_greater_degrevlex : public std::binary_function<ex, ex, bool> {
	const tinfo_t function_id;// = find_tinfo_key("function");
	const tinfo_t fderivative_id;// = find_tinfo_key("fderivative");
	const tinfo_t power_id;// = find_tinfo_key("power");
	const tinfo_t symbol_id;// = find_tinfo_key("symbol");
	const tinfo_t mul_id;// = find_tinfo_key("mul");
	const tinfo_t add_id;// = find_tinfo_key("add");
	const tinfo_t numeric_id;// = find_tinfo_key("numeric");
	const tinfo_t constant_id;// = find_tinfo_key("constant");
	const tinfo_t in_type_id;

	ex_is_greater_degrevlex(const tinfo_t type_id):
		function_id(find_tinfo_key("function")),
		fderivative_id(find_tinfo_key("fderivative")),
		power_id(find_tinfo_key("power")),
		symbol_id(find_tinfo_key("symbol")),
		mul_id(find_tinfo_key("mul")),
		add_id(find_tinfo_key("add")),
		numeric_id(find_tinfo_key("numeric")),
		constant_id(find_tinfo_key("constant")),
		in_type_id(type_id) {};

	bool operator() (const ex &lh, const ex &rh) const;
	int compare(const ex &lh, const ex &rh) const;
	int compare(const basic *lh, const basic *rh) const;
	// mul objects
	int compare_mul_symbol(const mul *lh, const symbol *rh) const;
	int compare_mul_power(const mul *lh, const power *rh) const;
	int compare_same_type_mul(const mul *lh, const mul *rh) const;
	// add objects
	int compare_add_symbol(const add *lh, const symbol *rh) const;
	int compare_add_power(const add *lh, const power *rh) const;
	int compare_add_mul(const add *lh, const mul *rh) const;
	int compare_same_type_add(const add *lh, const add *rh) const;
	// power objects
	int compare_power_symbol(const power *lh, const symbol *rh) const;
	int compare_same_type_power(const power *lh, const power *rh) const;
	// symbol objects
	int compare_same_type_symbol(const symbol *lh, const symbol *rh) const;
	// container objects
	template <template <class T, class = std::allocator<T> > class C>
	int compare_same_type_container(const container<C> *lh,const container<C> *rh) const;
	// function objects
	int compare_same_type_function(const function *lh, const function *rh) const;
	// fderivative objects
	int compare_same_type_fderivative(const fderivative *lh, const fderivative *rh) const;

	static const ex_is_greater_degrevlex& in_type(tinfo_t in_type_id) {
	  static ex_is_greater_degrevlex in_type[2] = {ex_is_greater_degrevlex(&add::tinfo_static),
						       ex_is_greater_degrevlex(&mul::tinfo_static)};
	        if (in_type_id == &mul::tinfo_static)
		        return in_type[1];
        	return in_type[0];
	}
};

// We have to define the following class to sort held expressions
// E.g. 3*x+2*x which does not get simplified to 5*x.
struct expair_is_greater_degrevlex : public std::binary_function<expair, expair, bool>
{
        const tinfo_t in_type_id;
        expair_is_greater_degrevlex(tinfo_t in_type):
	        in_type_id(in_type) {};
	bool operator() (const expair &lh, const expair &rh) const;

	static const expair_is_greater_degrevlex& in_type(tinfo_t in_type_id) {
        	static expair_is_greater_degrevlex in_type[2] = {expair_is_greater_degrevlex(&add::tinfo_static),expair_is_greater_degrevlex(&add::tinfo_static)};
		if (in_type_id == &mul::tinfo_static)
		        return in_type[1];
        	return in_type[0];
	}
};

struct expair_rest_is_greater_degrevlex : public std::binary_function<expair, expair, bool>
{
        const tinfo_t in_type_id;
        expair_rest_is_greater_degrevlex(tinfo_t in_type):
	        in_type_id(in_type) {};
	bool operator() (const expair &lh, const expair &rh) const {
	        return ex_is_greater_degrevlex::in_type(in_type_id)(lh.rest, rh.rest);
	}

	static const expair_rest_is_greater_degrevlex& in_type(tinfo_t in_type_id) {
        	static expair_rest_is_greater_degrevlex in_type[2] = {expair_rest_is_greater_degrevlex(&add::tinfo_static),
					expair_rest_is_greater_degrevlex(&mul::tinfo_static)};
		if (in_type_id == &mul::tinfo_static)
		        return in_type[1];
        	return in_type[0];
	}
};

} // namespace GiNaC

#endif // ndef __GINAC_ORDER_H__
