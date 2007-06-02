/*

frobenius.h

Copyright (C) 2007, David Harvey

Please see frobenius.cpp for licensing information.

*/


#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>


/*
Computes frobenius matrix for given p, to precision p^N, for the
hyperelliptic curve y^2 = Q(x), on the standard basis of cohomology.

PRECONDITIONS:
   p must be a prime > (2g+1)(2N-1).
   N >= 1.
   Degree of Q should be 2g+1 for some g >= 1.
   Q must be monic. The reduction of Q mod p must have no multiple roots.

RETURN VALUE:
   1 on success, in which case "output" holds the resulting 2g*2g matrix.
   0 if any of the above conditions are not satisfied (EXCEPTION: frobenius()
       will not check that p is prime. That's up to you.)

*/
int frobenius(NTL::mat_ZZ& output, const NTL::ZZ& p, int N, const NTL::ZZX& Q);


// ----------------------- end of file
