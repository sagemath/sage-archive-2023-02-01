// -*- c++ -*-
//*****************************************************************************
/** @file BooleEnv.cc
 *
 * @author Alexander Dreyer
 * @date 2006-03-06
 *
 * This file implements the class BooleEnv, which hold global structures.
 * a polynomial ring over Booleans.
 *
 * @par Copyright:
 *   (c) 2007 by The PolyBoRi Team
 *
 * @internal
 * @version \$Id$
 *
 * @par History:
 * @verbatim
 * $Log$
 * Revision 1.4  2009/07/23 19:41:06  dreyer
 * ADD: BooleRing::hash
 *
 * Revision 1.3  2008/04/29 09:01:52  bricken
 * + active_ring moved to BooleEnv.cc
 *
 * Revision 1.2  2007/12/13 16:18:07  dreyer
 * CHANGE: removed unnecessary friend declaration
 *
 * Revision 1.1  2007/12/13 15:53:49  dreyer
 * CHANGE: Ordering in BoolePolyRing again; BooleEnv manages active ring
 *
 * @endverbatim
**/
//*****************************************************************************


// load header file
# include "BooleEnv.h"
# include "BoolePolyRing.h"
#include "OrderedManager.h"

BEGIN_NAMESPACE_PBORI


//BooleEnv::ring_type active_ring;
// Note, direct access via friends, to  active_ring.pOrder and pMgt, because
// inline doesn't work up to now, because os the undefined type issue.

/// @todo needs inlining!!!
#ifdef PBORI_ENV_RING_NOTINLINED
BooleEnv::ring_type& BooleEnv::ring() {
  static BooleEnv::ring_type active_ring(1000, CTypes::lp, false);
  return active_ring;
}
#endif

BooleEnv::block_iterator
BooleEnv::blockBegin() {

  return ordering().blockBegin();
}

BooleEnv::block_iterator
BooleEnv::blockEnd() {

  return ordering().blockEnd();
}

void BooleEnv::appendBlock(idx_type idx) {

  ordering().appendBlock(idx);
}

void BooleEnv::clearBlocks() {

  ordering().clearBlocks();
}



BooleEnv::idx_type
BooleEnv::lastBlockStart() {
  return ring().lastBlockStart();
}




BooleEnv::manager_type& BooleEnv::manager() {
  return ring().manager(); }
BooleEnv::order_type& BooleEnv::ordering() {
  return  ring().ordering(); }




  /// Get empty decision diagram
BooleEnv::dd_type BooleEnv::zero() { return ring().zero(); }

  /// Get decision diagram with all variables negated
BooleEnv::dd_type BooleEnv::one() { return ring().one(); }

  /// Get number of ring variables the of active ring
BooleEnv::size_type BooleEnv::nVariables() {
  return manager().nVariables();
}





  /// Set name of variable with index idx
void
BooleEnv::setVariableName(idx_type idx, vartext_type varname) {
  ring().setVariableName(idx, varname);
}

  /// Get name of variable with index idx
BooleEnv::vartext_type
BooleEnv::getVariableName(idx_type idx){
  return ring().getVariableName(idx);
}


  /// Change order of current ring
void
BooleEnv::changeOrdering(ordercode_type code) {
    ring().changeOrdering(code);
}




  /// Get numerical code for current ordering
BooleEnv::ordercode_type BooleEnv::getOrderCode() {
  return ordering().getOrderCode();
}

  /// Get numerical code for current base ordering
  /// (the same for non-block orderings)
BooleEnv::ordercode_type BooleEnv::getBaseOrderCode() {
  return ordering().getBaseOrderCode();
}


void
BooleEnv::printInfo() {

  return ring().printInfo();
}



  /// Access idx-th variable of the active ring
BooleEnv::dd_type BooleEnv::variable(idx_type idx) {
  return manager().variable(idx);
}


  /// Access idx-th variable
BooleEnv::dd_type BooleEnv::persistentVariable(idx_type idx) {
    return manager().persistentVariable(idx);
  }


void BooleEnv::set(ring_type& theRing) { ring() = theRing; }



END_NAMESPACE_PBORI
