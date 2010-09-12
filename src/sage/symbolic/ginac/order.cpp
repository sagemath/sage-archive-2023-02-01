
#include "order.h"
#include "registrar.h"
#include "utils.h"

#include <cmath>
#include <iostream>

namespace GiNaC {


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
		} else {
			return lh->compare_same_type(*rh);
		}
	} else if (typeid_lh == fderivative_id) {
		//print fderivatives after everything
		return -1;
	} else if (typeid_lh == function_id) {
		//print functions before fderivatives, after anything else
		return typeid_rh == fderivative_id ? 1 : -1;
	} else if (typeid_lh == mul_id) {
		if (typeid_rh == power_id) {
			return compare_mul_power(
					static_cast<const mul*>(lh),
					static_cast<const power*>(rh));
		} else if (typeid_rh == symbol_id) {
			return compare_mul_symbol(
					static_cast<const mul*>(lh),
					static_cast<const symbol*>(rh));
		}
	} else if (typeid_lh == power_id) {
		if (typeid_rh == mul_id) {
			return -compare_mul_power(
					static_cast<const mul*>(rh),
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
		} else if (typeid_rh == power_id) {
			return -compare_power_symbol(
					static_cast<const power*>(rh),
					static_cast<const symbol*>(lh));
		}
	}
	//std::cout<<"comparing typeid's"<<std::endl;
	return (typeid_lh<typeid_rh ? 1 : -1);
}

int ex_is_greater_degrevlex::compare_mul_symbol(const mul *lh,
		const symbol *rh) const
{
	int cmpval;
	double tdeg;
	tdeg = lh->total_degree();
	if (tdeg == 1) {
		cmpval = compare(get_pointer(lh->seq[0].rest.bp), rh);
		if (cmpval != 0) {
			return cmpval;
		}
		cmpval = compare(_num1_p,
				get_pointer(lh->seq[0].coeff.bp));
		if (cmpval != 0) {
			return cmpval;
		}
		return -1;
	}
	return tdeg > 1 ? 1 : -1;
}

// compare this to a pow object
// first we compare degrees
// if equal we compare the first item in the sequence to the base in other
int ex_is_greater_degrevlex::compare_mul_power(const mul *lh,
		const power *rh) const
{
	double lh_deg = lh->total_degree();
	double rh_deg;
	numeric rh_exp;
	int cmpval = 0;
	if (is_a<numeric>(rh->exponent)) {
		rh_exp = ex_to<numeric>(rh->exponent);
		if (rh_exp.is_real()) {
			rh_deg = rh_exp.to_double();
		} else {
			rh_deg = std::sqrt(std::pow(rh_exp.real().to_double(), 2) + 
					std::pow(rh_exp.imag().to_double(), 2));
		}
		if (rh_deg != lh_deg)
			return lh_deg < rh_deg ? -1 : 1;
	} else {
		cmpval = compare(get_pointer(lh->seq[0].coeff.bp),
				get_pointer(rh->exponent.bp));
		if (cmpval != 0)
			return cmpval;
	}
	cmpval = compare(get_pointer(lh->seq[0].rest.bp),
			get_pointer(rh->basis.bp));
	if (cmpval != 0) {
		return cmpval;
	}
	if (lh->seq.size() == 1 && lh->overall_coeff.is_equal(_ex1))
		return 0;
	return 1;
}

int ex_is_greater_degrevlex::compare_same_type_mul(const mul *lh,
		const mul *rh) const
{
	int cmpval;

	// compare total degrees
	double lh_deg = lh->total_degree();
	double rh_deg = rh->total_degree();
	if (lh_deg != rh_deg)
		return lh_deg < rh_deg ? -1 : 1;

	// compare each item in lh to corresponding element in rh
	epvector::const_iterator cit1 = lh->seq.begin();
	epvector::const_iterator cit2 = rh->seq.begin();
	epvector::const_iterator last1 = lh->seq.end();
	epvector::const_iterator last2 = rh->seq.end();

	for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
		cmpval = ( cit1->compare(*cit2));
		if (cmpval != 0)
			return cmpval;
	}

	// compare sizes
	if (cit1 != last1) 
		return 1;
	else if (cit2 != last2)
		return -1;

	// compare overall_coeff
	cmpval = compare(get_pointer(lh->overall_coeff.bp),
			get_pointer(rh->overall_coeff.bp));
	if (cmpval!=0)
		return cmpval;

	return 0;
}

int ex_is_greater_degrevlex::compare_same_type_add(const add *lh,
		const add *rh) const
{
	int cmpval;

	// compare number of elements
	if (lh->seq.size() != rh->seq.size())
		return (lh->seq.size()<rh->seq.size()) ? -1 : 1;

	epvector::const_iterator cit1 = lh->seq.begin();
	epvector::const_iterator cit2 = rh->seq.begin();
	epvector::const_iterator last1 = lh->seq.end();
	epvector::const_iterator last2 = rh->seq.end();

	for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
		cmpval = compare(get_pointer(cit2->rest.bp),
				get_pointer(cit1->rest.bp));
		if (cmpval!=0) return cmpval;
	}

	GINAC_ASSERT(cit1==last1);
	GINAC_ASSERT(cit2==last2);

	// compare overall_coeff
	cmpval = compare(get_pointer(lh->overall_coeff.bp),
			get_pointer(rh->overall_coeff.bp));
	if (cmpval!=0)
		return cmpval;

	return 0;
}

int ex_is_greater_degrevlex::compare_power_symbol(const power *lh,
		const symbol *rh) const
{
	//std::cout<<"in compare_power_symbol"<<std::endl;
	int cmpval;
	cmpval = compare(get_pointer(lh->exponent.bp), _num1_p);
	if (cmpval != 0) {
		return cmpval;
	}
	return compare(get_pointer(lh->basis.bp), rh);
}

int ex_is_greater_degrevlex::compare_same_type_power(const power *lh,
		const power *rh) const
{
	int cmpval = compare(get_pointer(lh->exponent.bp),
				get_pointer(rh->exponent.bp));
	if (cmpval)
		return cmpval;
	else
		return compare(get_pointer(lh->basis.bp),
			get_pointer(rh->basis.bp));
}

int ex_is_greater_degrevlex::compare_same_type_symbol(const symbol *lh,
		const symbol *rh) const
{
	//std::cout<<"in compare_same symbol"<<std::endl;
	//std::cout<<"lh: ";
	//lh->dbgprint();
	//std::cout<<"rh: ";
	//rh->dbgprint();
	//std::cout<<std::endl;
	if (lh->serial==rh->serial) return 0;
	//std::cout<<"after if"<<std::endl;

// SAGE/Pynac: Sorting based on creation order doesn't work for Sage.
// instead we sort on variable name. -- William Stein

	//std::cout<<"comparing names: "<<(lh->name < rh->name)<<std::endl;
	return lh->name < (rh->name) ? 1 : -1;
}

} // namespace GiNaC
