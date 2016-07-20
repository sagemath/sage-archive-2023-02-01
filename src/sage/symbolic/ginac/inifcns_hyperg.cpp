/** @file inifcns_hyperg.cpp
 *
 *  Implementation of hypergeometric functions.
 *
 *  (C) 2016 Ralf Stephan <ralf@ark.in-berlin.de>
 */

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

#include "inifcns.h"
#include "ex.h"
#include "constant.h"
#include "infinity.h"
#include "numeric.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "symbol.h"
#include "pseries.h"
#include "utils.h"

#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <memory>

namespace GiNaC {


//////////
// Appell F1 function
//////////

static ex appellf1_evalf(const ex& a, const ex& b1, const ex& b2, const ex& c, const ex& x, const ex& y, PyObject* parent)
{
	/*if (is_exactly_a<numeric>(a) and is_exactly_a<numeric>(b1)
            and is_exactly_a<numeric>(b2) and is_exactly_a<numeric>(c)
            and is_exactly_a<numeric>(x) and is_exactly_a<numeric>(y)) {
                return appellf1(ex_to<numeric>(a), ex_to<numeric>(b1),
                                ex_to<numeric>(b2), ex_to<numeric>(c),
                                ex_to<numeric>(x), ex_to<numeric>(y));
	}*/
	return appell_F1(a, b1, b2, c, x, y).hold();
}

static ex appellf1_eval(const ex& a, const ex& b1, const ex& b2, const ex& c, const ex& x, const ex& y)
{
	if (is_exactly_a<numeric>(a) and is_exactly_a<numeric>(b1)
            and is_exactly_a<numeric>(b2) and is_exactly_a<numeric>(c)
            and is_exactly_a<numeric>(x) and is_exactly_a<numeric>(y)) {
                return appellf1_evalf(a, b1, b2, c, x, y, nullptr);
	}

	return appell_F1(a, b1, b2, c, x, y).hold();
}

static ex appellf1_deriv(const ex& a, const ex& b1, const ex& b2, const ex& c, const ex& x, const ex& y, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==4 || deriv_param==5);
	
	// d/dx F1(a,b1,b2,c,x,y) --> a*b1/c*F1(a+1,b1+1,b2,c+1,x,y)
	// d/dy F1(a,b1,b2,c,x,y) --> a*b2/c*F1(a+1,b1,b2+1,c+1,x,y)
        if (deriv_param == 4)
        	return mul(mul(mul(a, b1), pow(c, -1)), appell_F1(a+1, b1+1, b2, c+1, x, y));
        return mul(mul(mul(a, b2), pow(c, -1)), appell_F1(a+1, b1, b2+1, c+1, x, y));
}

REGISTER_FUNCTION(appell_F1, eval_func(appellf1_eval).
                        evalf_func(appellf1_evalf).
                        derivative_func(appellf1_deriv).
		        latex_name("\\operatorname{F_1}"));

} // namespace GiNaC
