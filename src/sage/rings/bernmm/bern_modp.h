/*
   bern_modp.cpp:  computing isolated Bernoulli numbers modulo p

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

#ifndef BERNMM_BERN_MODP_H
#define BERNMM_BERN_MODP_H


namespace bernmm {


/*
   Returns B_k mod p, in [0, p), or -1 if B_k is not p-integral.

   PRECONDITIONS:
      2 <= p < NTL_SP_BOUND, p prime
      k >= 0
*/
long bern_modp(long p, long k);


/*
   Exported for testing.
*/
long _bern_modp_powg(long p, double pinv, long k);
long _bern_modp_pow2(long p, double pinv, long k);


};


#endif

// end of file ================================================================
