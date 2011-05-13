// -*- c++ -*-
//*****************************************************************************
/** @file COrderedIter.h
 *
 * @author Alexander Dreyer
 * @date 2006-09-06
 *
 * This file defines an iterator, which respects the current ordering.
 *
 * @par Copyright:
 *   (c) 2006 by The PolyBoRi Team
**/
//*****************************************************************************

#ifndef COrderedIter_h_
#define COrderedIter_h_

// include basic definitions
#include "pbori_defs.h"
#include "pbori_algo.h"


#include "BoolePolynomial.h"
  //#include "OrderedManager.h"
#include "CDelayedTermIter.h"
#include "CBidirectTermIter.h"
#include <algorithm>

#include "CStackSelector.h"
#include "CTermGenerator.h"


BEGIN_NAMESPACE_PBORI


template <class NavigatorType>
class CAbstractStackBase {
public:
  typedef NavigatorType navigator;

  typedef CAbstractStackBase<NavigatorType> self;
  typedef CTermStackBase<NavigatorType, self> iterator_core;
  typedef boost::shared_ptr<iterator_core> core_pointer;

  virtual void increment() = 0;
  virtual core_pointer copy() const = 0;

  virtual ~CAbstractStackBase() {}
};



template <class StackType>
class CWrappedStack:
  public StackType {
public:
  typedef StackType base;
  typedef CWrappedStack<StackType> self;

  typedef typename base::navigator navigator;

  typedef typename base::iterator_core iterator_core;
  typedef boost::shared_ptr<iterator_core> core_pointer;

  template <class MgrType>
  CWrappedStack(navigator navi, const MgrType& mgr):
    base(navi, mgr) {
    base::init();
  }
  CWrappedStack(): base() {}
  CWrappedStack(const self& rhs): base(rhs) {}


  core_pointer copy() const {
    return core_pointer(new self(*this));
  }

};


// template<class SequenceType>
// void get_term(BooleMonomial& monom, const SequenceType& seq) {

//   typename SequenceType::const_reverse_iterator start(seq.rbegin()),
//     finish(seq.rend());

//   while (start != finish){
//     monom.changeAssign(*start);
//     ++start;
//   }
// }


// template<class SequenceType>
// void get_term(BooleExponent& termexp, const SequenceType& seq) {

//   termexp.reserve(seq.deg());
//   typename SequenceType::const_iterator start(seq.begin()),
//     finish(seq.end());

//   while (start != finish){
//     termexp.push_back(*start);
//     ++start;
//   }
// }


// template<class SequenceType>
// void get_term(typename CTypes::size_type& termdeg, const SequenceType& seq) {

//   termdeg = seq.deg();
// }

template <class NavigatorType, class MonomType>
class COrderedIter:
  public boost::iterator_facade<
  COrderedIter<NavigatorType, MonomType>,
  MonomType, std::forward_iterator_tag, MonomType
  > {

public:

  typedef COrderedIter<NavigatorType, MonomType> self;
  typedef CAbstractStackBase<NavigatorType> stack_base;
  typedef CTermStackBase<NavigatorType, stack_base> iterator_core;

  /// Type for functional, which generates actual term, for current path
  typedef CTermGenerator<MonomType> term_generator;

  typedef typename iterator_core::const_iterator const_iterator;
  typedef typename iterator_core::const_reverse_iterator
  const_reverse_iterator;
  typedef typename iterator_core::size_type size_type;
  typedef typename iterator_core::deg_type deg_type;
  typedef typename iterator_core::idx_type idx_type;


  /// Fix type of direct iterator
  typedef NavigatorType navigator;

  // Store shared pointer of iterator
  typedef boost::shared_ptr<iterator_core> core_pointer;

  /// Extract plain Boolean type
  typedef bool bool_type;

  // Constructor
  COrderedIter(core_pointer rhs,
               const term_generator & getTerm):
    m_getTerm(getTerm), p_iter(rhs) {}

  // Destructor
  ~COrderedIter() {}

  bool equal(const self& rhs) const {
    return  p_iter->equal(*rhs.p_iter); }

  /// Incrementation
  void increment() {
    if (!p_iter.unique()) {
      core_pointer tmp(p_iter->copy());
      p_iter = tmp;
    }

    p_iter->increment();
  }

  /// Determine whether term is one (without explicit constructing)
  bool_type isOne() const { return p_iter->isOne(); }

  /// Determine whether term is zero (without explicit constructing)
  bool_type isZero() const { return p_iter->isZero(); }

  /// Check, whether end of iteration is reached
  bool_type isEnd() const { return isZero(); }

  /// Dereferencing operation
  MonomType dereference() const {

    return m_getTerm(*p_iter);
  }

  const_iterator begin() const { return p_iter->begin(); }
  const_iterator end() const { return p_iter->end(); }
  const_reverse_iterator rbegin() const { return p_iter->rbegin(); }
  const_reverse_iterator rend() const { return p_iter->rend(); }

  deg_type deg() const { return p_iter->deg(); }
  idx_type firstIndex() const { return *begin(); }

  /// Get navigator of term start
  navigator navigation() const {
    return p_iter->navigation();
  }

protected:
  /// The functional which defines the dereferecing operation
  term_generator m_getTerm;

  /// A shared pointer to the stack, which carries the current path
  core_pointer p_iter;
};


template <class OrderType, class NavigatorType, class MonomType>
class CGenericOrderedIter:
  public COrderedIter<NavigatorType, MonomType> {
public:
  typedef CAbstractStackBase<NavigatorType> stack_base;
  typedef typename CStackSelector<OrderType, NavigatorType, stack_base>::type
  ordered_iter_base;
  typedef CWrappedStack<ordered_iter_base> ordered_iter_type;

  typedef COrderedIter<NavigatorType, MonomType> base;
  typedef typename base::iterator_core iterator_core;
  typedef typename base::core_pointer core_pointer;

  typedef typename base::term_generator term_generator;

  template <class MgrType>
  CGenericOrderedIter(NavigatorType navi, const MgrType& gen):
    base( core_pointer(new ordered_iter_type(navi, gen) ), gen) {}
  CGenericOrderedIter(): base( core_pointer(new ordered_iter_type()),
                               term_generator() ) {}

  CGenericOrderedIter(const CGenericOrderedIter& rhs): base(rhs) {}
};

template <class OrderType, class NavigatorType>
class CGenericOrderedIter<OrderType, NavigatorType, BooleExponent> :
  public COrderedIter<NavigatorType, BooleExponent> {
public:
  typedef CAbstractStackBase<NavigatorType> stack_base;
  typedef typename CStackSelector<OrderType, NavigatorType, stack_base>::type
  ordered_iter_base;
  typedef CWrappedStack<ordered_iter_base> ordered_iter_type;

  typedef COrderedIter<NavigatorType, BooleExponent> base;
  typedef typename base::iterator_core iterator_core;
  typedef typename base::core_pointer core_pointer;

  typedef typename base::term_generator term_generator;

  template <class MgrType>
  CGenericOrderedIter(NavigatorType navi, const MgrType& mgr):
    base( core_pointer(new ordered_iter_type(navi, mgr)),
                       term_generator() ) {}

  CGenericOrderedIter(): base( core_pointer(new ordered_iter_type()),
                              term_generator() ) {}

  CGenericOrderedIter(const CGenericOrderedIter& rhs): base(rhs) {}
};

END_NAMESPACE_PBORI

#endif
