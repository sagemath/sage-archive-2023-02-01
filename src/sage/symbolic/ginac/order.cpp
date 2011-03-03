/** @file order.cpp
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

#include "order.h"
#include "registrar.h"
#include "utils.h"

#include <cmath>
#include <iostream>

namespace GiNaC {

/** What SAGE does for printing:
 To print multivariate polynomials, SAGE uses "reversed" (bigger terms first) 'degrevlex' order
 i.e. "reversed" graded reversed lexicographic order, the lexicographic order is defined below.
 The variables are ordered according to their "creation" order
 i.e PR.<x,a,t> = PolynomialRing(QQ) gives x > a > t
 and x+t+a is printed x + a + t
*/

/** Returns true if lhex > rhex for 'degrevlex' of SAGE
 i.e. graded reversed lexicographic order
*/
bool ex_is_greater_degrevlex::operator() (const ex &lhex, const ex &rhex) const {
        //std::cout<<"in ex_is_greater_degrevlex::operator()"<<std::endl;
	//std::cout<<"lh: ";
	//lhex.dbgprint();
	//std::cout<<std::endl<<"rh: ";
	//rhex.dbgprint();
	//std::cout<<std::endl;
	const basic *lh = get_pointer(lhex.bp);
	const basic *rh = get_pointer(rhex.bp);
	return compare(lh, rh) == 1;
}

int ex_is_greater_degrevlex::compare (const ex &lhex, const ex &rhex) const {
	const basic *lh = get_pointer(lhex.bp);
	const basic *rh = get_pointer(rhex.bp);
	return compare(lh, rh);
}

bool expair_is_greater_degrevlex::operator() (const expair &lhex, const expair &rhex) const
{
	const basic *lh = get_pointer(lhex.rest.bp);
	const basic *rh = get_pointer(rhex.rest.bp);

	// compare rests
        int cmpval = ex_is_greater_degrevlex::in_type(in_type_id).compare(lh, rh);
	if (cmpval != 0) {
	        return cmpval == 1;
	}

	// compare coeffs which are numerics
	//std::cout << "comparing coeffs" << std::endl;
	//lhex.coeff.bp->dbgprinttree();
	//rhex.coeff.bp->dbgprinttree();
	numeric lh_coeff = ex_to<numeric>(lhex.coeff);
	numeric rh_coeff = ex_to<numeric>(rhex.coeff);
	// strange...
	//std::cout<<lh->compare_same_type(rh)<<std::endl;
	//return lh->compare_same_type(rh) == 1;
	double lh_deg = 1;
	double rh_deg = 1;
		if (lh_coeff.is_real()) {
			lh_deg = lh_coeff.to_double();
		} else {
			lh_deg = std::sqrt(std::pow(lh_coeff.real().to_double(), 2) + 
					std::pow(lh_coeff.imag().to_double(), 2));
		}
		if (rh_coeff.is_real()) {
			rh_deg = rh_coeff.to_double();
		} else {
			rh_deg = std::sqrt(std::pow(rh_coeff.real().to_double(), 2) + 
					std::pow(rh_coeff.imag().to_double(), 2));
		}

	return lh_deg > rh_deg;	
}

/** Comparison functions:
 They should implement 'degrevlex' of SAGE
 i.e. graded reversed lexicographic order
 The lexicographic order should depend on the "creation" order of variables ?
 Or should it be a "natural" one on strings : a > b > ... > x > y > ... ?
*/

/** Return values for comparison functions :
 compare(a,b) should return :
 -1 if a < b
 0 if a == b
 1 if a > b
 as <=> in Perl and GiNaC internal functions
*/

int ex_is_greater_degrevlex::compare(const basic *lh, const basic *rh) const {
	const tinfo_t typeid_lh = lh->tinfo();
	const tinfo_t typeid_rh = rh->tinfo();

	if (typeid_rh==typeid_lh) {
		if (typeid_rh == mul_id) {
			return compare_same_type_mul(
					static_cast<const mul*>(lh),
					static_cast<const mul*>(rh));
		} else if (typeid_rh == add_id) {
			return compare_same_type_add(
					static_cast<const add*>(lh),
					static_cast<const add*>(rh));
		} else if (typeid_rh == symbol_id) {
			return compare_same_type_symbol(
					static_cast<const symbol*>(lh),
					static_cast<const symbol*>(rh));
		} else if (typeid_rh == power_id) {
			return compare_same_type_power(
					static_cast<const power*>(lh),
					static_cast<const power*>(rh));
		} else if (typeid_rh == function_id) {
			return compare_same_type_function(
					static_cast<const function*>(lh),
					static_cast<const function*>(rh));
		} else if (typeid_rh == fderivative_id) {
			return compare_same_type_fderivative(
					static_cast<const fderivative*>(lh),
					static_cast<const fderivative*>(rh));
		} else {
		        // using GiNaC functions by default
			return lh->compare_same_type(*rh);
		}
	// at present numerics are combined into overall_coefficient
	// even when hold parameter is used
	} else if (typeid_lh == numeric_id) {
	 	//print numerics after anything else
	        return -1;
	} else if (typeid_rh == numeric_id) {
	 	//print numerics after anything else
	        return 1;
	} else if (typeid_lh == constant_id) {
	 	//print constants after anything else
	        return -1;
	} else if (typeid_rh == constant_id) {
	 	//print constants after anything else
	        return 1;
	} else if (typeid_lh == fderivative_id) {
		//print fderivatives after everything else
		return -1;
	} else if (typeid_rh == fderivative_id) {
		//print fderivatives after everything esle
		return 1;
	} else if (typeid_lh == function_id) {
		//print functions before fderivatives, after anything else
	        return -1;
	} else if (typeid_rh == function_id) {
		//print functions before fderivatives, after anything else
	        return 1;
	} else if (typeid_lh == mul_id) {
		if (typeid_rh == power_id) {
			return compare_mul_power(
					static_cast<const mul*>(lh),
					static_cast<const power*>(rh));
		} else if (typeid_rh == symbol_id) {
			return compare_mul_symbol(
					static_cast<const mul*>(lh),
					static_cast<const symbol*>(rh));
		} else if (typeid_rh == add_id) {
			return -compare_add_mul(
					static_cast<const add*>(rh),
					static_cast<const mul*>(lh));
		}
	} else if (typeid_lh == add_id) {
		if (typeid_rh == power_id) {
			return compare_add_power(
					static_cast<const add*>(lh),
					static_cast<const power*>(rh));
		} else if (typeid_rh == symbol_id) {
			return compare_add_symbol(
					static_cast<const add*>(lh),
					static_cast<const symbol*>(rh));
		} else if (typeid_rh == mul_id) {
			return compare_add_mul(
					static_cast<const add*>(lh),
					static_cast<const mul*>(rh));
		}
	} else if (typeid_lh == power_id) {
		if (typeid_rh == mul_id) {
			return -compare_mul_power(
					static_cast<const mul*>(rh),
					static_cast<const power*>(lh));
		} else if (typeid_rh == add_id) {
			return -compare_add_power(
					static_cast<const add*>(rh),
					static_cast<const power*>(lh));
		} else if (typeid_rh == symbol_id) {
			return compare_power_symbol(
					static_cast<const power*>(lh),
					static_cast<const symbol*>(rh));
		}
	} else if (typeid_lh == symbol_id) {
		if (typeid_rh == mul_id) {
			return -compare_mul_symbol(
					static_cast<const mul*>(rh),
					static_cast<const symbol*>(lh));
		} else if (typeid_rh == add_id) {
			return -compare_add_symbol(
					static_cast<const add*>(rh),
					static_cast<const symbol*>(lh));
		} else if (typeid_rh == power_id) {
			return -compare_power_symbol(
					static_cast<const power*>(rh),
					static_cast<const symbol*>(lh));
		}
        }
	std::cout<<"comparing typeid's"<<std::endl;
	return (typeid_lh<typeid_rh ? -1 : 1);
}

// compare a mul and a symbol objects
// same behavior within mul and add objects:
// the total degree of the symbol is 1
// check total degree of mul then compare smallest item to symbol
int ex_is_greater_degrevlex::compare_mul_symbol(const mul *lh,
		const symbol *rh) const
{
  //std::cout<<"comparing mul and symbol"<<std::endl;
  //lh->dbgprint();
  //rh->dbgprint();
	int cmpval;

	double tdeg;
	tdeg = lh->total_degree();

	if (tdeg != 1)
	        return tdeg > 1 ? 1 : -1;
	
	const expair smallest_item = lh->get_sorted_seq()->back();

	// compare bases
	cmpval = compare(get_pointer(smallest_item.rest.bp), rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare exponents
	cmpval = -compare(get_pointer(smallest_item.coeff.bp),_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh->seq.size() == 1 && lh->overall_coeff.is_equal(_ex1))
		return 0;

	// there is something of total degree 0 in front of the mul object
	return 1;
}


// compare a mul and a pow objects
// same behavior wihtin mul and add objects:
// first we compare total degrees
// if equal we compare the smallest basis in the sequence to the basis in other
// then their exponents
int ex_is_greater_degrevlex::compare_mul_power(const mul *lh,
		const power *rh) const
{
	int cmpval;

	double lh_deg = lh->total_degree();
	double rh_deg = 1;
	numeric rh_exp;
	if (is_a<numeric>(rh->exponent)) {
		rh_exp = ex_to<numeric>(rh->exponent);
		if (rh_exp.is_real()) {
			rh_deg = rh_exp.to_double();
		} else {
			rh_deg = std::sqrt(std::pow(rh_exp.real().to_double(), 2) + 
					std::pow(rh_exp.imag().to_double(), 2));
		}
	}
	if (rh_deg != lh_deg)
		return lh_deg < rh_deg ? -1 : 1;
	// same total degree

	// smallest item is at the end of the sorted sequence
	const expair smallest_item = lh->get_sorted_seq()->back();
	
	// compare bases
	cmpval = compare(get_pointer(smallest_item.rest.bp),
			get_pointer(rh->basis.bp));
	if (cmpval != 0) {
		return cmpval;
	}

	// compare exponents
	cmpval = -compare(get_pointer(smallest_item.coeff.bp),
		        get_pointer(rh->exponent.bp));
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh->seq.size() == 1 && lh->overall_coeff.is_equal(_ex1))
		return 0;

	// there is something of total degree 0 in front of the mul object
	return 1;
}

// compare two mul objects
// same behavior wihtin mul and add objects:
// first we compare total degrees
// if equal we compare the basis of the smallest items
// then their exponents
// and so on
int ex_is_greater_degrevlex::compare_same_type_mul(const mul *lh,
		const mul *rh) const
{
  //std::cout<<"comparing mul and mul"<<std::endl;
  //lh->dbgprint();
  //rh->dbgprint();
	int cmpval;
	
	// compare total degrees
	double lh_deg = lh->total_degree();
	double rh_deg = rh->total_degree();
        //std::cout<<"degree "<<lh_deg<<" and "<<rh_deg<<std::endl;
	if (lh_deg != rh_deg)
		return lh_deg < rh_deg ? -1 : 1;
	
	// compare each item in lh to corresponding element in rh
	const epvector *sorted_seq1 = lh->get_sorted_seq();
	const epvector *sorted_seq2 = rh->get_sorted_seq();
	epvector::const_reverse_iterator cit1 = sorted_seq1->rbegin();
	epvector::const_reverse_iterator cit2 = sorted_seq2->rbegin();
	epvector::const_reverse_iterator last1 = sorted_seq1->rend();
	epvector::const_reverse_iterator last2 = sorted_seq2->rend();

	for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
		// compare bases
		cmpval = compare(get_pointer(cit1->rest.bp),get_pointer(cit2->rest.bp));
		if (cmpval != 0) {
			return cmpval;
		}

		// compare exponents
	        cmpval = -compare(get_pointer(cit1->coeff.bp),get_pointer(cit2->coeff.bp));
		if (cmpval != 0) {
			return cmpval;
		}
	}

	// compare sizes
	if (cit1 != last1) 
		return 1;
	else if (cit2 != last2)
		return -1;

	// compare overall_coeff
	cmpval = compare(get_pointer(lh->overall_coeff.bp),
			get_pointer(rh->overall_coeff.bp));

	return cmpval;
}

// compare an add and a symbol objects
// same behavior wihtin mul and add objects:
// the coefficient of the symbol is 1
int ex_is_greater_degrevlex::compare_add_symbol(const add *lh,
		const symbol *rh) const
{
  //std::cout<<"comparing add and symbol"<<std::endl;
  //lh->dbgprint();
  //rh->dbgprint();
	int cmpval;

	const expair biggest_item = lh->get_sorted_seq()->front();

	// compare bases
	cmpval = ex_is_greater_degrevlex::in_type(&add::tinfo_static).compare(get_pointer(biggest_item.rest.bp), rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare coefficients
	cmpval = compare(get_pointer(biggest_item.coeff.bp),_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh->seq.size() == 1 && lh->overall_coeff.is_equal(_ex0))
		return 0;

	// there is something at the end of the add object
	return 1;
}

// compare an add and a mul objects
// same behavior within mul and add objects:
// the coefficient of the mul object is 1
int ex_is_greater_degrevlex::compare_add_mul(const add *lh,
		const mul *rh) const
{
        int cmpval;
	const expair biggest_item = lh->get_sorted_seq()->front();

	// compare bases
	cmpval = ex_is_greater_degrevlex::in_type(&add::tinfo_static).compare(get_pointer(biggest_item.rest.bp), rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare coefficients
	cmpval = compare(get_pointer(biggest_item.coeff.bp),_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh->seq.size() == 1 && lh->overall_coeff.is_equal(_ex0))
		return 0;

	// there is something at the end of the object
	return 1;
}

// compare an add and a pow objects
// same behavior wihtin mul and add objects:
// the coefficient of the power object is 1
int ex_is_greater_degrevlex::compare_add_power(const add *lh,
		const power *rh) const
{
        int cmpval;
	const expair biggest_item = lh->get_sorted_seq()->front();

	// compare bases
	cmpval = ex_is_greater_degrevlex::in_type(&add::tinfo_static).compare(get_pointer(biggest_item.rest.bp), rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare coefficients
	cmpval = compare(get_pointer(biggest_item.coeff.bp),_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh->seq.size() == 1 && lh->overall_coeff.is_equal(_ex0))
		return 0;

	// there is something at the end of the object
	return 1;
}

// compare two add objects
// same behavior wihtin mul and add objects:
// first we compare the basis of the biggest items
// then their coefficients
// and so on
int ex_is_greater_degrevlex::compare_same_type_add(const add *lh,
		const add *rh) const
{
	int cmpval;

	const epvector *sorted_seq1 = lh->get_sorted_seq();
	const epvector *sorted_seq2 = rh->get_sorted_seq();
	epvector::const_iterator cit1 = sorted_seq1->begin();
	epvector::const_iterator cit2 = sorted_seq2->begin();
	epvector::const_iterator last1 = sorted_seq1->end();
	epvector::const_iterator last2 = sorted_seq2->end();

	for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
		// compare bases
		cmpval = ex_is_greater_degrevlex::in_type(&add::tinfo_static).compare(get_pointer(cit1->rest.bp),get_pointer(cit2->rest.bp));
		if (cmpval != 0) {
			return cmpval;
		}

		// compare coefficients
	        cmpval = compare(get_pointer(cit1->coeff.bp),get_pointer(cit2->coeff.bp));
		if (cmpval != 0) {
			return cmpval;
		}
	}

	// compare sizes
	if (cit1 != last1) 
		return 1;
	else if (cit2 != last2)
		return -1;

	// compare overall_coeff
	cmpval = compare(get_pointer(lh->overall_coeff.bp),
			get_pointer(rh->overall_coeff.bp));
	
	return cmpval;
}

// compare a power and a symbol objects
// does not behave the same way inside add or mul
// in mul object:
// first we compare bases
// then exponents
// in add object:
// first exponents
// then bases
int ex_is_greater_degrevlex::compare_power_symbol(const power *lh,
		const symbol *rh) const
{
        int cmpval;

	//std::cout<<"in compare_power_symbol"<<std::endl;
        if (in_type_id == &mul::tinfo_static) {
	        cmpval = compare(get_pointer(lh->basis.bp), rh);
		if (cmpval != 0)
        		  return cmpval;
		
		cmpval = compare(get_pointer(lh->exponent.bp), _num1_p);
		
		return cmpval;
	}

	// We are in an add object
	if (is_a<numeric>(lh->exponent)) {
		numeric lh_exp = ex_to<numeric>(lh->exponent);
		double lh_deg;
		if (lh_exp.is_real()) {
			lh_deg = lh_exp.to_double();
		} else {
			lh_deg = std::sqrt(std::pow(lh_exp.real().to_double(), 2) + 
					std::pow(lh_exp.imag().to_double(), 2));
		}
		if (lh_deg != 1)
			return lh_deg < 1 ? -1 : 1;
	}

	cmpval = compare(get_pointer(lh->basis.bp), rh);

	return cmpval;
}

// compare two power objects
// does not behave the same way inside add or mul
// in mul object:
// fist we compare bases
// then exponents
// in add object:
// first exponents
// then bases
int ex_is_greater_degrevlex::compare_same_type_power(const power *lh,
	        const power *rh) const
{
        int cmpval;
        if (in_type_id == &mul::tinfo_static) {
	        cmpval = compare(get_pointer(lh->basis.bp), get_pointer(rh->basis.bp));
		if (cmpval != 0)
        		  return cmpval;
		cmpval = compare(get_pointer(lh->exponent.bp), get_pointer(rh->exponent.bp));
		
		return cmpval;
	}

	double lh_deg = 1;
	double rh_deg = 1;
	if (is_a<numeric>(lh->exponent)) {
		numeric lh_exp = ex_to<numeric>(lh->exponent);
		if (lh_exp.is_real()) {
			lh_deg = lh_exp.to_double();
		} else {
			lh_deg = std::sqrt(std::pow(lh_exp.real().to_double(), 2) + 
					std::pow(lh_exp.imag().to_double(), 2));
		}
	}
	if (is_a<numeric>(rh->exponent)) {
		numeric rh_exp = ex_to<numeric>(rh->exponent);
		if (rh_exp.is_real()) {
			rh_deg = rh_exp.to_double();
		} else {
			rh_deg = std::sqrt(std::pow(rh_exp.real().to_double(), 2) + 
					std::pow(rh_exp.imag().to_double(), 2));
		}
	}
	if (rh_deg != lh_deg)
		return lh_deg < rh_deg ? -1 : 1;

	cmpval = compare(get_pointer(lh->basis.bp),
			 get_pointer(rh->basis.bp));
	if (cmpval != 0)
        	return cmpval;

	if (is_a<numeric>(lh->exponent) && is_a<numeric>(rh->exponent))
	        return 0;
	return compare(get_pointer(lh->exponent.bp), get_pointer(rh->exponent.bp));
}

// compare two symbol objects
// same behavior wihtin mul and add objects:
// we compare names
int ex_is_greater_degrevlex::compare_same_type_symbol(const symbol *lh,
		const symbol *rh) const
{
        //std::cout<<"in compare_same symbol"<<std::endl;
	//std::cout<<"lh: ";
	//lh->dbgprint();
	//std::cout<<"rh: ";
	//rh->dbgprint();
	//std::cout<<std::endl;

        // Weird because Sage sorts based on creation order.
        // SAGE/Pynac: Sorting based on creation order doesn't work for Sage.
        // instead we sort on variable name. -- William Stein

	//std::cout<<"comparing names: "<<(lh->name < rh->name)<<std::endl;
	/* Reversed ordering on char encoding (i.e. "a" < "b") then length (i.e. "abc" < "abcd")
	   i.e. "x" > "y" and "xyz" > "xyzt" */
	if (lh->serial==rh->serial) return 0;
	//std::cout<<"after if"<<std::endl;
	return lh->name < rh->name ? 1 : -1;

	//return lh->serial < rh->serial ? 1 : -1;
	//return -lh->compare_same_type(*rh);
}

// compare two containers of the same type
template <template <class T, class = std::allocator<T> > class C>
int ex_is_greater_degrevlex::compare_same_type_container(const container<C> *lh,
							 const container<C> *rh) const
{
        typename C<ex>::const_iterator it1 = lh->seq.begin(), it1end = lh->seq.end(),
	                      it2 = rh->seq.begin(), it2end = rh->seq.end();

	while (it1 != it1end && it2 != it2end) {
	        int cmpval = compare(get_pointer(it1->bp), get_pointer(it2->bp));
		if (cmpval)
			return cmpval;
		++it1; ++it2;
	}

	return (it1 == it1end) ? (it2 == it2end ? 0 : -1) : 1;
}

// compare two function objects
// same behavior wihtin mul and add objects:
// we compare names
int ex_is_greater_degrevlex::compare_same_type_function(const function *lh,
		const function *rh) const
{

	//std::cout<<"comparing names: "<<(lh->get_name() < rh->get_name())<<std::endl;
	/* Reversed ordering on char encoding (i.e. "a" < "b") then length (i.e. "abc" < "abcd")
	   i.e. "x" > "y" and "xyz" > "xyzt" */
        if (lh->serial==rh->serial) //return 0;
	        return compare_same_type_container(lh,rh);
	//std::cout<<"after if"<<std::endl;
	return lh->get_name() < rh->get_name() ? 1 : -1;	

	//return lh->serial < rh->serial ? 1 : -1;
	//return -lh->compare_same_type(*rh);
}

// compare two fderivative objects
// same behavior wihtin mul and add objects:
// we compare names
int ex_is_greater_degrevlex::compare_same_type_fderivative(const fderivative *lh,
		const fderivative *rh) const
{
        int cmpval = compare_same_type_function(lh, rh);
	if (cmpval != 0)
	        return cmpval;
	if (lh->parameter_set != rh->parameter_set)
	        return lh->parameter_set < rh->parameter_set ? -1 : 1;
	return 0;
}

} // namespace GiNaC
