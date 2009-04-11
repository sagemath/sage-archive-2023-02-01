/** @file function.h
 *
 *  Interface to class of symbolic functions. */

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

#ifndef __GINAC_FUNCTION_H__
#define __GINAC_FUNCTION_H__

#include <string>
#include <vector>

// CINT needs <algorithm> to work properly with <vector>
#include <algorithm>

#include "exprseq.h"

// the following lines have been generated for max. 14 parameters
#define DECLARE_FUNCTION_1P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 1; \
template<typename T1> const GiNaC::function NAME(const T1 & p1) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1)); \
}

#define DECLARE_FUNCTION_2P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 2; \
template<typename T1, typename T2> const GiNaC::function NAME(const T1 & p1, const T2 & p2) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2)); \
}

#define DECLARE_FUNCTION_3P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 3; \
template<typename T1, typename T2, typename T3> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3)); \
}

#define DECLARE_FUNCTION_4P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 4; \
template<typename T1, typename T2, typename T3, typename T4> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4)); \
}

#define DECLARE_FUNCTION_5P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 5; \
template<typename T1, typename T2, typename T3, typename T4, typename T5> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5)); \
}

#define DECLARE_FUNCTION_6P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 6; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6)); \
}

#define DECLARE_FUNCTION_7P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 7; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7)); \
}

#define DECLARE_FUNCTION_8P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 8; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7, const T8 & p8) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7), GiNaC::ex(p8)); \
}

#define DECLARE_FUNCTION_9P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 9; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7, const T8 & p8, const T9 & p9) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7), GiNaC::ex(p8), GiNaC::ex(p9)); \
}

#define DECLARE_FUNCTION_10P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 10; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7, const T8 & p8, const T9 & p9, const T10 & p10) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7), GiNaC::ex(p8), GiNaC::ex(p9), GiNaC::ex(p10)); \
}

#define DECLARE_FUNCTION_11P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 11; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7, const T8 & p8, const T9 & p9, const T10 & p10, const T11 & p11) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7), GiNaC::ex(p8), GiNaC::ex(p9), GiNaC::ex(p10), GiNaC::ex(p11)); \
}

#define DECLARE_FUNCTION_12P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 12; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7, const T8 & p8, const T9 & p9, const T10 & p10, const T11 & p11, const T12 & p12) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7), GiNaC::ex(p8), GiNaC::ex(p9), GiNaC::ex(p10), GiNaC::ex(p11), GiNaC::ex(p12)); \
}

#define DECLARE_FUNCTION_13P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 13; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7, const T8 & p8, const T9 & p9, const T10 & p10, const T11 & p11, const T12 & p12, const T13 & p13) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7), GiNaC::ex(p8), GiNaC::ex(p9), GiNaC::ex(p10), GiNaC::ex(p11), GiNaC::ex(p12), GiNaC::ex(p13)); \
}

#define DECLARE_FUNCTION_14P(NAME) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = 14; \
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14> const GiNaC::function NAME(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6, const T7 & p7, const T8 & p8, const T9 & p9, const T10 & p10, const T11 & p11, const T12 & p12, const T13 & p13, const T14 & p14) { \
	return GiNaC::function(NAME##_SERIAL::serial, GiNaC::ex(p1), GiNaC::ex(p2), GiNaC::ex(p3), GiNaC::ex(p4), GiNaC::ex(p5), GiNaC::ex(p6), GiNaC::ex(p7), GiNaC::ex(p8), GiNaC::ex(p9), GiNaC::ex(p10), GiNaC::ex(p11), GiNaC::ex(p12), GiNaC::ex(p13), GiNaC::ex(p14)); \
}


// end of generated lines

#define REGISTER_FUNCTION(NAME,OPT) \
unsigned NAME##_SERIAL::serial = \
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

// the following lines have been generated for max. 14 parameters
typedef ex (* eval_funcp_1)(const ex &);
typedef ex (* eval_funcp_2)(const ex &, const ex &);
typedef ex (* eval_funcp_3)(const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_4)(const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* eval_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);

typedef ex (* evalf_funcp_1)(const ex &);
typedef ex (* evalf_funcp_2)(const ex &, const ex &);
typedef ex (* evalf_funcp_3)(const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_4)(const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* evalf_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);

typedef ex (* conjugate_funcp_1)(const ex &);
typedef ex (* conjugate_funcp_2)(const ex &, const ex &);
typedef ex (* conjugate_funcp_3)(const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_4)(const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* conjugate_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);

typedef ex (* real_part_funcp_1)(const ex &);
typedef ex (* real_part_funcp_2)(const ex &, const ex &);
typedef ex (* real_part_funcp_3)(const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_4)(const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* real_part_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);

typedef ex (* imag_part_funcp_1)(const ex &);
typedef ex (* imag_part_funcp_2)(const ex &, const ex &);
typedef ex (* imag_part_funcp_3)(const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_4)(const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* imag_part_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);

typedef ex (* derivative_funcp_1)(const ex &, unsigned);
typedef ex (* derivative_funcp_2)(const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_3)(const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_4)(const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);
typedef ex (* derivative_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, unsigned);

typedef ex (* power_funcp_1)(const ex &, const ex &);
typedef ex (* power_funcp_2)(const ex &, const ex &, const ex &);
typedef ex (* power_funcp_3)(const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_4)(const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);
typedef ex (* power_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &);

typedef ex (* series_funcp_1)(const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_2)(const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_3)(const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_4)(const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);
typedef ex (* series_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const relational &, int, unsigned);

typedef void (* print_funcp_1)(const ex &, const print_context &);
typedef void (* print_funcp_2)(const ex &, const ex &, const print_context &);
typedef void (* print_funcp_3)(const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_4)(const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_5)(const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_6)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_7)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_8)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_9)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_10)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_11)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_12)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_13)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);
typedef void (* print_funcp_14)(const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const ex &, const print_context &);

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
// the following lines have been generated for max. 14 parameters
    function_options & eval_func(eval_funcp_1 e);
    function_options & eval_func(eval_funcp_2 e);
    function_options & eval_func(eval_funcp_3 e);
    function_options & eval_func(eval_funcp_4 e);
    function_options & eval_func(eval_funcp_5 e);
    function_options & eval_func(eval_funcp_6 e);
    function_options & eval_func(eval_funcp_7 e);
    function_options & eval_func(eval_funcp_8 e);
    function_options & eval_func(eval_funcp_9 e);
    function_options & eval_func(eval_funcp_10 e);
    function_options & eval_func(eval_funcp_11 e);
    function_options & eval_func(eval_funcp_12 e);
    function_options & eval_func(eval_funcp_13 e);
    function_options & eval_func(eval_funcp_14 e);

    function_options & evalf_func(evalf_funcp_1 ef);
    function_options & evalf_func(evalf_funcp_2 ef);
    function_options & evalf_func(evalf_funcp_3 ef);
    function_options & evalf_func(evalf_funcp_4 ef);
    function_options & evalf_func(evalf_funcp_5 ef);
    function_options & evalf_func(evalf_funcp_6 ef);
    function_options & evalf_func(evalf_funcp_7 ef);
    function_options & evalf_func(evalf_funcp_8 ef);
    function_options & evalf_func(evalf_funcp_9 ef);
    function_options & evalf_func(evalf_funcp_10 ef);
    function_options & evalf_func(evalf_funcp_11 ef);
    function_options & evalf_func(evalf_funcp_12 ef);
    function_options & evalf_func(evalf_funcp_13 ef);
    function_options & evalf_func(evalf_funcp_14 ef);

    function_options & conjugate_func(conjugate_funcp_1 d);
    function_options & conjugate_func(conjugate_funcp_2 d);
    function_options & conjugate_func(conjugate_funcp_3 d);
    function_options & conjugate_func(conjugate_funcp_4 d);
    function_options & conjugate_func(conjugate_funcp_5 d);
    function_options & conjugate_func(conjugate_funcp_6 d);
    function_options & conjugate_func(conjugate_funcp_7 d);
    function_options & conjugate_func(conjugate_funcp_8 d);
    function_options & conjugate_func(conjugate_funcp_9 d);
    function_options & conjugate_func(conjugate_funcp_10 d);
    function_options & conjugate_func(conjugate_funcp_11 d);
    function_options & conjugate_func(conjugate_funcp_12 d);
    function_options & conjugate_func(conjugate_funcp_13 d);
    function_options & conjugate_func(conjugate_funcp_14 d);

    function_options & real_part_func(real_part_funcp_1 d);
    function_options & real_part_func(real_part_funcp_2 d);
    function_options & real_part_func(real_part_funcp_3 d);
    function_options & real_part_func(real_part_funcp_4 d);
    function_options & real_part_func(real_part_funcp_5 d);
    function_options & real_part_func(real_part_funcp_6 d);
    function_options & real_part_func(real_part_funcp_7 d);
    function_options & real_part_func(real_part_funcp_8 d);
    function_options & real_part_func(real_part_funcp_9 d);
    function_options & real_part_func(real_part_funcp_10 d);
    function_options & real_part_func(real_part_funcp_11 d);
    function_options & real_part_func(real_part_funcp_12 d);
    function_options & real_part_func(real_part_funcp_13 d);
    function_options & real_part_func(real_part_funcp_14 d);

    function_options & imag_part_func(imag_part_funcp_1 d);
    function_options & imag_part_func(imag_part_funcp_2 d);
    function_options & imag_part_func(imag_part_funcp_3 d);
    function_options & imag_part_func(imag_part_funcp_4 d);
    function_options & imag_part_func(imag_part_funcp_5 d);
    function_options & imag_part_func(imag_part_funcp_6 d);
    function_options & imag_part_func(imag_part_funcp_7 d);
    function_options & imag_part_func(imag_part_funcp_8 d);
    function_options & imag_part_func(imag_part_funcp_9 d);
    function_options & imag_part_func(imag_part_funcp_10 d);
    function_options & imag_part_func(imag_part_funcp_11 d);
    function_options & imag_part_func(imag_part_funcp_12 d);
    function_options & imag_part_func(imag_part_funcp_13 d);
    function_options & imag_part_func(imag_part_funcp_14 d);

    function_options & derivative_func(derivative_funcp_1 d);
    function_options & derivative_func(derivative_funcp_2 d);
    function_options & derivative_func(derivative_funcp_3 d);
    function_options & derivative_func(derivative_funcp_4 d);
    function_options & derivative_func(derivative_funcp_5 d);
    function_options & derivative_func(derivative_funcp_6 d);
    function_options & derivative_func(derivative_funcp_7 d);
    function_options & derivative_func(derivative_funcp_8 d);
    function_options & derivative_func(derivative_funcp_9 d);
    function_options & derivative_func(derivative_funcp_10 d);
    function_options & derivative_func(derivative_funcp_11 d);
    function_options & derivative_func(derivative_funcp_12 d);
    function_options & derivative_func(derivative_funcp_13 d);
    function_options & derivative_func(derivative_funcp_14 d);

    function_options & power_func(power_funcp_1 d);
    function_options & power_func(power_funcp_2 d);
    function_options & power_func(power_funcp_3 d);
    function_options & power_func(power_funcp_4 d);
    function_options & power_func(power_funcp_5 d);
    function_options & power_func(power_funcp_6 d);
    function_options & power_func(power_funcp_7 d);
    function_options & power_func(power_funcp_8 d);
    function_options & power_func(power_funcp_9 d);
    function_options & power_func(power_funcp_10 d);
    function_options & power_func(power_funcp_11 d);
    function_options & power_func(power_funcp_12 d);
    function_options & power_func(power_funcp_13 d);
    function_options & power_func(power_funcp_14 d);

    function_options & series_func(series_funcp_1 s);
    function_options & series_func(series_funcp_2 s);
    function_options & series_func(series_funcp_3 s);
    function_options & series_func(series_funcp_4 s);
    function_options & series_func(series_funcp_5 s);
    function_options & series_func(series_funcp_6 s);
    function_options & series_func(series_funcp_7 s);
    function_options & series_func(series_funcp_8 s);
    function_options & series_func(series_funcp_9 s);
    function_options & series_func(series_funcp_10 s);
    function_options & series_func(series_funcp_11 s);
    function_options & series_func(series_funcp_12 s);
    function_options & series_func(series_funcp_13 s);
    function_options & series_func(series_funcp_14 s);

    template <class Ctx> function_options & print_func(print_funcp_1 p)
    {
    	test_and_set_nparams(1);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_2 p)
    {
    	test_and_set_nparams(2);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_3 p)
    {
    	test_and_set_nparams(3);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_4 p)
    {
    	test_and_set_nparams(4);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_5 p)
    {
    	test_and_set_nparams(5);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_6 p)
    {
    	test_and_set_nparams(6);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_7 p)
    {
    	test_and_set_nparams(7);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_8 p)
    {
    	test_and_set_nparams(8);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_9 p)
    {
    	test_and_set_nparams(9);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_10 p)
    {
    	test_and_set_nparams(10);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_11 p)
    {
    	test_and_set_nparams(11);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_12 p)
    {
    	test_and_set_nparams(12);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_13 p)
    {
    	test_and_set_nparams(13);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }
    template <class Ctx> function_options & print_func(print_funcp_14 p)
    {
    	test_and_set_nparams(14);
    	set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp(p));
    	return *this;
    }

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

	// python function calls
	function_options & eval_func(PyObject* e);
	function_options & evalf_func(PyObject* e);
	function_options & conjugate_func(PyObject* e);
	function_options & real_part_func(PyObject* e);
	function_options & imag_part_func(PyObject* e);
	function_options & derivative_func(PyObject* e);
	function_options & power_func(PyObject* e);
	function_options & series_func(PyObject* e);

	function_options & set_return_type(unsigned rt, tinfo_t rtt=NULL);
	function_options & do_not_evalf_params();
	function_options & remember(unsigned size, unsigned assoc_size=0,
	                            unsigned strategy=remember_strategies::delete_never);
	function_options & overloaded(unsigned o);
	function_options & set_symmetry(const symmetry & s);

	std::string get_name() const { return name; }
	unsigned get_nparams() const { return nparams; }

	void set_python_func() { python_func = true; }

	void set_print_latex_func(PyObject* f);
	void set_print_dflt_func(PyObject* f);
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

	bool python_func;

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
	// the following lines have been generated for max. 14 parameters
    function(unsigned ser, const ex & param1);
    function(unsigned ser, const ex & param1, const ex & param2);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12, const ex & param13);
    function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12, const ex & param13, const ex & param14);

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
	int compare(const basic &other) const;
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
	bool lookup_remember_table(ex & result) const;
	void store_remember_table(ex const & result) const;
public:
	static std::vector<function_options> & registered_functions();
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

