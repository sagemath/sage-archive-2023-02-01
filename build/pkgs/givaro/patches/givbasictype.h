#ifndef _GIVARO_BASICTYPE_H_
#define _GIVARO_BASICTYPE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givbasictype.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givbasictype.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
#include "givaro/givconfig.h"

#include <stdlib.h> // for size_t
#ifdef MACOSX
#  include <sys/types.h> // needed on MacOS X 10.5 for uint type
#endif

// -- Neutral type: definition of zero and one
class Neutral {
public:
  static Neutral zero;
  static Neutral one;
  inline operator int() const { return _val; }
  inline int operator==( const Neutral& n) const { return _val==n._val; }
  inline int operator!=( const Neutral& n) const { return _val!=n._val; }
private:
  Neutral( int val ) : _val(val) {};
  int _val;
};

// -- Used to build no initialized object as static object
class givNoInit {};
// -- Used to call cstor without copy
class givNoCopy {};
// -- Used to call cstor with copy
class givWithCopy {};

#endif
