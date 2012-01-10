/** @file normal.h
 *
 *  This file defines several functions that work on univariate and
 *  multivariate polynomials and rational functions.
 *  These functions include polynomial quotient and remainder, GCD and LCM
 *  computation, square-free factorization and rational function normalization. */

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

#ifndef __GINAC_NORMAL_H__
#define __GINAC_NORMAL_H__

#include "lst.h"

namespace GiNaC {

class ex;
class symbol;

// Quotient q(x) of polynomials a(x) and b(x) in Q[x], so that a(x)=b(x)*q(x)+r(x)
extern ex quo(const ex &a, const ex &b, const ex &x, bool check_args = true);

// Remainder r(x) of polynomials a(x) and b(x) in Q[x], so that a(x)=b(x)*q(x)+r(x)
extern ex rem(const ex &a, const ex &b, const ex &x, bool check_args = true);

// Decompose rational function a(x)=N(x)/D(x) into Q(x)+R(x)/D(x) with degree(R, x) < degree(D, x)
extern ex decomp_rational(const ex &a, const ex &x);

// Pseudo-remainder of polynomials a(x) and b(x) in Q[x]
extern ex prem(const ex &a, const ex &b, const ex &x, bool check_args = true);

// Pseudo-remainder of polynomials a(x) and b(x) in Q[x]
extern ex sprem(const ex &a, const ex &b, const ex &x, bool check_args = true);

// Exact polynomial division of a(X) by b(X) in Q[X] (quotient returned in q), returns false when exact division fails
extern bool divide(const ex &a, const ex &b, ex &q, bool check_args = true);

// Polynomial GCD in Z[X], cofactors are returned in ca and cb, if desired
extern ex gcd(const ex &a, const ex &b, ex *ca = NULL, ex *cb = NULL, bool check_args = true);

// Polynomial LCM in Z[X]
extern ex lcm(const ex &a, const ex &b, bool check_args = true);

// Square-free factorization of a polynomial a(x)
extern ex sqrfree(const ex &a, const lst &l = lst());

// Square-free partial fraction decomposition of a rational function a(x)
extern ex sqrfree_parfrac(const ex & a, const symbol & x);

// Collect common factors in sums.
extern ex collect_common_factors(const ex & e);

// Resultant of two polynomials e1,e2 with respect to symbol s.
extern ex resultant(const ex & e1, const ex & e2, const ex & s);

} // namespace GiNaC

#endif // ndef __GINAC_NORMAL_H__
