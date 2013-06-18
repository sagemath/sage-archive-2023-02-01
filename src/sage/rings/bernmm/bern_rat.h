/*
   bern_rat.h:  multi-modular algorithm for computing Bernoulli numbers

   Copyright (C) 2008, 2009, David Harvey

   This file is part of the bernmm package (version 1.1).

   bernmm is released under a BSD-style license. See the README file in
   the source distribution for details.
*/

#ifndef BERNMM_BERN_RAT_H
#define BERNMM_BERN_RAT_H

#include <gmp.h>


namespace bernmm {


/*
   Returns B_k as a rational number, stored at _res_.

   k must be >= 0.

   Uses _num_threads_ threads. (If USE_THREADS is not #defined, the code just
   uses one thread.)
*/
void bern_rat(mpq_t res, long k, int num_threads);



};    // end namespace


#endif

// end of file ================================================================
