#include "ex.h"
#include "basic.h"
#include "mul.h"
#include "power.h"
#include "add.h"
#include "symbol.h"

namespace GiNaC {
	
	struct ex_is_greater_degrevlex : public std::binary_function<ex, ex, bool> {
		const tinfo_t function_id;// = find_tinfo_key("function");
		const tinfo_t fderivative_id;// = find_tinfo_key("fderivative");
		const tinfo_t power_id;// = find_tinfo_key("power");
		const tinfo_t symbol_id;// = find_tinfo_key("symbol");
		const tinfo_t mul_id;// = find_tinfo_key("mul");
		const tinfo_t add_id;// = find_tinfo_key("add");

		ex_is_greater_degrevlex():
			function_id(find_tinfo_key("function")),
			fderivative_id(find_tinfo_key("fderivative")),
			power_id(find_tinfo_key("power")),
			symbol_id(find_tinfo_key("symbol")),
			mul_id(find_tinfo_key("mul")),
			add_id(find_tinfo_key("add")) {};
		bool operator() (const ex &lh, const ex &rh) const;
		int compare(const basic *lh, const basic *rh) const;
		int compare_mul_symbol(const mul *lh, const symbol *rh) const;
		int compare_mul_power(const mul *lh, const power *rh) const;
		int compare_same_type_mul(const mul *lh, const mul *rh) const;
		int compare_same_type_add(const add *lh, const add *rh) const;
		int compare_power_symbol(const power *lh, const symbol *rh) const;
		int compare_same_type_power(const power *lh, const power *rh) const;
		int compare_same_type_symbol(const symbol *lh, const symbol *rh) const;
	};

struct expair_rest_is_greater_degrevlex : public std::binary_function<expair, expair, bool> {
	bool operator()(const expair &lh, const expair &rh) const { return ex_is_greater_degrevlex()(lh.rest, rh.rest); }
};

}
