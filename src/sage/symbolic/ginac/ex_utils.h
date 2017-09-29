#ifndef __PYNAC_EXUTILS_H__
#define __PYNAC_EXUTILS_H___

#include "ex.h"
#include "numeric.h"

namespace GiNaC {


// utility functions

// wrapper functions around member functions
inline size_t nops(const ex & thisex)
{ return thisex.nops(); }

inline ex expand(const ex & thisex, unsigned options = 0)
{ return thisex.expand(options); }

inline ex conjugate(const ex & thisex)
{ return thisex.conjugate(); }

inline ex real_part(const ex & thisex)
{ return thisex.real_part(); }

inline ex imag_part(const ex & thisex)
{ return thisex.imag_part(); }

inline bool has(const ex & thisex, const ex & pattern, unsigned options = 0)
{ return thisex.has(pattern, options); }

inline bool find(const ex & thisex, const ex & pattern, lst & found)
{ return thisex.find(pattern, found); }

inline bool is_polynomial(const ex & thisex, const ex & vars)
{ return thisex.is_polynomial(vars); }

inline numeric degree(const ex & thisex, const ex & s)
{ return thisex.degree(s); }

inline numeric ldegree(const ex & thisex, const ex & s)
{ return thisex.ldegree(s); }

inline ex coeff(const ex & thisex, const ex & s, int n=1)
{ return thisex.coeff(s, n); }

inline ex numer(const ex & thisex)
{ return thisex.numer(); }

inline ex denom(const ex & thisex)
{ return thisex.denom(); }

inline ex numer_denom(const ex & thisex)
{ return thisex.numer_denom(); }

inline ex normal(const ex & thisex, int level=0, bool noexpand_combined=false,
                bool noexpand_numer=true)
{ return thisex.normal(level, noexpand_combined, noexpand_numer); }

inline ex to_rational(const ex & thisex, lst & repl_lst)
{ return thisex.to_rational(repl_lst); }

inline ex to_rational(const ex & thisex, exmap & repl)
{ return thisex.to_rational(repl); }

inline ex to_polynomial(const ex & thisex, exmap & repl)
{ return thisex.to_polynomial(repl); }

inline ex to_polynomial(const ex & thisex, lst & repl_lst)
{ return thisex.to_polynomial(repl_lst); }

inline ex collect(const ex & thisex, const ex & s, bool distributed = false)
{ return thisex.collect(s, distributed); }

inline ex eval(const ex & thisex, int level = 0)
{ return thisex.eval(level); }

inline ex evalf(const ex & thisex, int level = 0, PyObject* parent=nullptr)
{ return thisex.evalf(level, parent); }

inline ex diff(const ex & thisex, const symbol & s, unsigned nth = 1)
{ return thisex.diff(s, nth); }

inline ex series(const ex & thisex, const ex & r, int order, unsigned options = 0)
{ return thisex.series(r, order, options); }

inline bool match(const ex & thisex, const ex & pattern, lst & repl_lst)
{ return thisex.match(pattern, repl_lst); }

inline ex symmetrize(const ex & thisex)
{ return thisex.symmetrize(); }

inline ex symmetrize(const ex & thisex, const lst & l)
{ return thisex.symmetrize(l); }

inline ex antisymmetrize(const ex & thisex)
{ return thisex.antisymmetrize(); }

inline ex antisymmetrize(const ex & thisex, const lst & l)
{ return thisex.antisymmetrize(l); }

inline ex symmetrize_cyclic(const ex & thisex)
{ return thisex.symmetrize_cyclic(); }

inline ex symmetrize_cyclic(const ex & thisex, const lst & l)
{ return thisex.symmetrize_cyclic(l); }

inline ex op(const ex & thisex, size_t i)
{ return thisex.op(i); }

inline ex lhs(const ex & thisex)
{ return thisex.lhs(); }

inline ex rhs(const ex & thisex)
{ return thisex.rhs(); }

inline bool is_zero(const ex & thisex)
{ return thisex.is_zero(); }

inline void swap(ex & e1, ex & e2)
{ e1.swap(e2); }

inline ex subs(const ex & thisex, const exmap & m, unsigned options = 0)
{ return thisex.subs(m, options); }

inline ex subs(const ex & thisex, const lst & ls, const lst & lr, unsigned options = 0)
{ return thisex.subs(ls, lr, options); }

inline ex subs(const ex & thisex, const ex & e, unsigned options = 0)
{ return thisex.subs(e, options); }


} // namespace GiNaC


// Specializations of Standard Library algorithms
namespace std {

/** Specialization of std::swap() for ex objects. */
template <>
inline void swap(GiNaC::ex &a, GiNaC::ex &b)
{
	a.swap(b);
}

/** Specialization of std::iter_swap() for vector<ex> iterators. */
template <>
inline void iter_swap(vector<GiNaC::ex>::iterator i1, vector<GiNaC::ex>::iterator i2)
{
	i1->swap(*i2);
}

/** Specialization of std::iter_swap() for list<ex> iterators. */
template <>
inline void iter_swap(list<GiNaC::ex>::iterator i1, list<GiNaC::ex>::iterator i2)
{
	i1->swap(*i2);
}

} // namespace std

#endif // ndef __PYNAC_EXUTILS_H___

