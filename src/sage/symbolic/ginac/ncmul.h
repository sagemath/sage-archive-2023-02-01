/** @file ncmul.h
 *
 *  Interface to GiNaC's non-commutative products of expressions. */

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

#ifndef __GINAC_NCMUL_H__
#define __GINAC_NCMUL_H__

#include "exprseq.h"

namespace GiNaC {

/** Non-commutative product of expressions. */
class ncmul : public exprseq
{
	GINAC_DECLARE_REGISTERED_CLASS(ncmul, exprseq)

	friend class power;
	friend ex reeval_ncmul(const exvector & v);
	friend ex hold_ncmul(const exvector & v);

// member functions

	// other constructors
public:
	ncmul(const ex & lh, const ex & rh);
	ncmul(const ex & f1, const ex & f2, const ex & f3);
	ncmul(const ex & f1, const ex & f2, const ex & f3,
	      const ex & f4);
	ncmul(const ex & f1, const ex & f2, const ex & f3,
	      const ex & f4, const ex & f5);
	ncmul(const ex & f1, const ex & f2, const ex & f3,
	      const ex & f4, const ex & f5, const ex & f6);
	ncmul(const exvector & v, bool discardable=false);
	ncmul(std::auto_ptr<exvector> vp);

	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const {return 50;}
	bool info(unsigned inf) const;
	int degree(const ex & s) const;
	int ldegree(const ex & s) const;
	ex expand(unsigned options=0) const;
	ex coeff(const ex & s, int n=1) const;
	ex eval(int level=0) const;
	ex evalm() const;
	exvector get_free_indices() const;
	ex thiscontainer(const exvector & v) const;
	ex thiscontainer(std::auto_ptr<exvector> vp) const;
	ex conjugate() const;
	ex real_part() const;
	ex imag_part() const;

protected:
	ex derivative(const symbol & s) const;
	unsigned return_type() const;
	tinfo_t return_type_tinfo() const;
	
	// new virtual functions which can be overridden by derived classes
	// none

	// non-virtual functions in this class
protected:
	void do_print(const print_context & c, unsigned level) const;
	void do_print_csrc(const print_context & c, unsigned level) const;
	size_t count_factors(const ex & e) const;
	void append_factors(exvector & v, const ex & e) const;
	std::auto_ptr<exvector> expandchildren(unsigned options) const;
public:
	const exvector & get_factors() const;
};

// friend funtions 

ex reeval_ncmul(const exvector & v);
ex hold_ncmul(const exvector & v);

} // namespace GiNaC

#endif // ndef __GINAC_NCMUL_H__
