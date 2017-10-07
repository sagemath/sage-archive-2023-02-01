/*
  This is coxgroup.cpp
  
  Coxeter version 3.0 Copyright (C) 2009 Mike Hansen
  See file main.cpp for full copyright notice
*/

#include "sage.h"

#include "error.h"

namespace sage {

  void interval(List<CoxWord>& list, coxgroup::CoxGroup& W, const CoxWord& g, const CoxWord& h)

  /* 
     Returns a list of the elements in the Bruhat interval between g and h.
     Note that this assumes that g and h are in order.
   */
  {
    if (not W.inOrder(g,h)) {
      return;
    }

    W.extendContext(h);
    
    CoxNbr x = W.contextNumber(g);
    CoxNbr y = W.contextNumber(h);
    
    BitMap b(W.contextSize());
    W.extractClosure(b,y);

    BitMap::ReverseIterator b_rend = b.rend();
    List<CoxNbr> res(0);

    for (BitMap::ReverseIterator i = b.rbegin(); i != b_rend; ++i)
      if (not W.inOrder(x,*i)) {
        BitMap bi(W.contextSize());
        W.extractClosure(bi,*i);
        CoxNbr z = *i; // andnot will invalidate iterator
        b.andnot(bi);
        b.setBit(z);   // otherwise the decrement will not be correct
      } else
        res.append(*i);

    schubert::NFCompare nfc(W.schubert(),W.ordering());
    Permutation a(res.size());
    sortI(res,nfc,a);
    
    list.setSize(0);
    for (size_t j = 0; j < res.size(); ++j) {
      CoxWord w(0);
      W.schubert().append(w, res[a[j]]);
      list.append(w);
    }
   
    return;
  }

}
