/* 
 * File:   assume.cpp
 * Author: Ralf Stephan <ralf@ark.in-berlin.de>
 * 
 * Created on August 5, 2015, 9:31 AM
 */

#include "ex.h"
#include "operators.h"
#include "symbol.h"
#include "function.h"
#include "relational.h"
#include "utils.h"
#include "flags.h"
#include "infoflagbase.h"

namespace GiNaC {

void assume(ex rel) {
        // It was already checked that rel is a relational.
        relational r = ex_to<relational>(rel);
        if (r.the_operator() == relational::equal
                or r.the_operator() == relational::not_equal)
                return;
        ex df = (r.lhs() - r.rhs()).expand();
        if (r.the_operator() == relational::greater)
                df.set_domain(domain::positive);
}

void assume(ex x, char* flag_desc) {
        if (!strcmp(flag_desc, "integer"))
                x.set_domain(domain::integer);
        else if (!strcmp(flag_desc, "real"))
                x.set_domain(domain::real);
        else if (!strcmp(flag_desc, "complex"))
                x.set_domain(domain::complex);
}

void forget(ex rel) {
        // It was already checked that rel is a relational.
        relational r = ex_to<relational>(rel);
        if (r.the_operator() == relational::equal
                or r.the_operator() == relational::not_equal)
                return;
        ex df = (r.lhs() - r.rhs()).expand();
        df.set_domain(domain::complex);
}

void forget(ex x, char* flag_desc) {
        x.set_domain(domain::complex);
}

} // namespace GiNaC
