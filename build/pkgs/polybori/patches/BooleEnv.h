// -*- c++ -*-
//*****************************************************************************
/** @file BooleEnv.h
 *
 * @author Alexander Dreyer
 * @date 2006-03-06
 *
 * This file the class BooleEnv, where handles global (static) strucutres of
 * PolyBoRi.
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
 * Revision 1.2  2008/04/29 09:02:10  bricken
 * + active_ring moved to BooleEnv.cc
 *
 * Revision 1.1  2007/12/13 15:53:48  dreyer
 * CHANGE: Ordering in BoolePolyRing again; BooleEnv manages active ring
 *
 * @endverbatim
**/
//*****************************************************************************


// load PolyBoRi settings
# include "pbori_defs.h"
#include <vector>

#ifndef BooleEnv_h_
#define BooleEnv_h_

// Obey stricter dependence of Sun Studio compiler
// todo: resolve dependency
#if 1 // defined(__SUNPRO_C) || defined(__SUNPRO_CC)
# define PBORI_ENV_RING_NOTINLINED
#endif

BEGIN_NAMESPACE_PBORI


class BoolePolyRing;
class CDynamicOrderBase;

/** @class BooleEnv
 * @brief This class is just a wrapper accessing global structures.
 *
 *
 **/


class BooleEnv:
  public CTypes::orderenums_type, public CTypes::compenums_type,
  public CTypes::auxtypes_type {

 public:
  //-------------------------------------------------------------------------
  // types definitions
  //-------------------------------------------------------------------------

  /// generic access to current type
  typedef BooleEnv self;

  /// generic access to base type
  typedef CTypes::orderenums_type base;

  /// @name adopt global type definitions
  //@{
  typedef CTypes::ordercode_type ordercode_type;
  typedef CTypes::manager_type manager_type;
  typedef CTypes::manager_reference manager_reference;
  typedef CTypes::manager_ptr manager_ptr;
  typedef CTypes::dd_type dd_type;
  typedef CTypes::vartext_type vartext_type;
  //@}

  /// Type for block indices
  typedef std::vector<idx_type> block_idx_type;

  /// Type for block iterators
  typedef block_idx_type::const_iterator block_iterator;

  //-------------------------------------------------------------------------
  // constructors and destructor
  //-------------------------------------------------------------------------

  /// Explicitely mention ordercodes' enumeration
  using base::ordercodes;

  /// Access idx-th variable of the active ring
  static dd_type variable(idx_type idx);

  /// Access idx-th variable
  static dd_type persistentVariable(idx_type idx);

  /// Get numerical code for current ordering
  static ordercode_type getOrderCode();

  /// Get numerical code for current base ordering
  /// (the same for non-block orderings)
  static ordercode_type getBaseOrderCode();

  /// Get empty decision diagram
  static dd_type zero();

  /// Get decision diagram with all variables negated
  static dd_type one();

  /// Get number of ring variables the of active ring
  static size_type nVariables();

  typedef BoolePolyRing ring_type;

  typedef CDynamicOrderBase order_type;
#ifdef PBORI_ENV_RING_NOTINLINED
  static ring_type& ring();
#else
  static ring_type& ring() {
    static BooleEnv::ring_type active_ring(1000, CTypes::lp, false);
    return active_ring;
  }
#endif
  static manager_type& manager();
  static order_type& ordering();
  /// Set name of variable with index idx
  static void setVariableName(idx_type idx, vartext_type varname);

  /// Get name of variable with index idx
  static vartext_type getVariableName(idx_type idx);

  /// @name interface for block orderings
  //@{
  static block_iterator blockBegin();
  static block_iterator blockEnd();
  static void appendBlock(idx_type idx);
  static void clearBlocks();

  static idx_type lastBlockStart();
  //@}

  /// Change order of current ring
  static void changeOrdering(ordercode_type code);

  static void printInfo();

  static void set(ring_type& theRing);


protected:


};

///please use BooleEnv::ring()

END_NAMESPACE_PBORI

#endif // of #ifndef BooleEnv_h_
