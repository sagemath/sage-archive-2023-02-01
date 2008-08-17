/*
   bern_rat.h:  multi-modular algorithm for computing Bernoulli numbers

   Copyright (C) 2008, David Harvey

   This file is part of the bernmm package (version 1.0).

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
