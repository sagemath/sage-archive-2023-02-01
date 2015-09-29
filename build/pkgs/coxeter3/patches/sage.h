/*
  This is interactive.h
  
  Coxeter version 3.0 Copyright (C) 2009 Mike Hansen
  See file main.cpp for full copyright notice
*/

#ifndef SAGE_H /* guard against multiple inclusions */
#define SAGE_H

#include "globals.h"
#include "coxgroup.h"
#include "coxtypes.h"
#include "schubert.h"

namespace sage {
  using namespace globals;
}

/******** function declarations **********************************************/

namespace sage {
  void interval(List<CoxWord>& result, coxgroup::CoxGroup& W, const CoxWord& g, const CoxWord& h);
}

#endif
