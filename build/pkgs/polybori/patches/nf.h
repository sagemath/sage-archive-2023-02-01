/*
 *  nf.h
 *  PolyBoRi
 *
 *  Created by Michael Brickenstein on 25.04.06.
 *  Copyright 2006 The PolyBoRi Team. See LICENSE file.
 *
 */
#include <algorithm>
#include <vector>
#include <queue>
#include "groebner_alg.h"
#if HAVE_M4RI
extern "C"{
#include "m4ri/m4ri.h"

#ifndef __M4RI_TWOPOW
#define __M4RI_TWOPOW TWOPOW
#endif
}
#endif
#ifndef PBORI_GB_NF_H
#define PBORI_GB_NF_H
BEGIN_NAMESPACE_PBORIGB


void drawmatrix(mzd_t* mat, const char* filename);

Polynomial add_up_exponents(const std::vector<Exponent>& vec);
Polynomial add_up_monomials(const std::vector<Monomial>& res_vec);


int select_no_deg_growth(const ReductionStrategy& strat, const Monomial& m);




class LMLessCompare{
public:
  bool operator() (const Polynomial& p1, const Polynomial& p2){
    return p1.lead()<p2.lead();
  }
};

const int SLIMGB_SIMPLEST=0;
template<int variant> class SlimgbReduction{
private:
  GroebnerStrategy* strat;
  std::priority_queue<Polynomial, std::vector<Polynomial>, LMLessCompare> to_reduce;
  public:
  std::vector<Polynomial> result;

  SlimgbReduction(GroebnerStrategy& strat){
    this->strat=&strat;
  }
  SlimgbReduction(){}
  void addPolynomial(const Polynomial& p);
  void reduce();
  //return zero at the end
  Polynomial nextResult();
};
template <int variant> void SlimgbReduction<variant>::addPolynomial(const Polynomial& p){
  if (!(p.isZero())){
    to_reduce.push(p);
  }
}
template <int variant> Polynomial SlimgbReduction<variant>::nextResult(){
  if (result.size()==0) return Polynomial();
  Polynomial res=result.back();
  result.pop_back();
  return res;
}
typedef SlimgbReduction<SLIMGB_SIMPLEST> slimgb_reduction_type;
std::vector<Polynomial> parallel_reduce(std::vector<Polynomial> inp, GroebnerStrategy& strat, int average_steps, double delay_f);
Polynomial red_tail(const ReductionStrategy& strat, Polynomial p);
Polynomial red_tail_short(const ReductionStrategy& strat, Polynomial p);
Polynomial nf3(const ReductionStrategy& strat, Polynomial p, Monomial rest_lead);
Polynomial nf3_short(const ReductionStrategy& strat, Polynomial p);
Polynomial ll_red_nf(const Polynomial& p,const BooleSet& reductors);

Polynomial ll_red_nf_noredsb(const Polynomial& p,const BooleSet& reductors);
Polynomial add_up_polynomials(const std::vector<Polynomial>& vec);
Polynomial plug_1(const Polynomial& p, const MonomialSet& m_plus_ones);
MonomialSet mod_mon_set(const MonomialSet& as, const MonomialSet &vs);
std::vector<Polynomial> gauss_on_polys(const std::vector<Polynomial>& orig_system);
Polynomial ll_red_nf_noredsb_single_recursive_call(const Polynomial& p,const BooleSet& reductors);
Polynomial cheap_reductions(const ReductionStrategy& strat, Polynomial p);
END_NAMESPACE_PBORIGB
#endif
