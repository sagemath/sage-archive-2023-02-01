#ifndef _GALREP_INCLUDE_
#define _GALREP_INCLUDE_

/*
    Copyright 2009 Andrew V. Sutherland

    This file is part of galrep.

    galrep is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    galrep is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with smalljac.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "gmp.h"

// NB: The int datatype is assumed to hold at least 32-bits

#define GALREP_MAX_ELL				59						// maximum prime ell
#define GALREP_ELL_COUNT				17						// # primes <= GALREP_MAX_ELL
#define GALREP_MAX_CCS				200						// upper bound on # of conjugacy classes spanning determinants in GL(2,Z/ellZ) for ell <= MAX_ELL
#define GALREP_MAX_GENERATORS			12

#define GALREP_BAD_ELL					-1						// ell out of range or not prime
#define GALREP_BAD_ERRBND				-2						// current errbnd=100 is fixed
#define GALREP_FILE_NOT_FOUND			-3						// couldn't find specified file (or default file)
#define GALREP_FILE_READ_ERROR			-4						// i/o error or unexpected eof
#define GALREP_BAD_ECDATA				-5						// invalid data found in elliptic curve data file
#define GALREP_BAD_GL2DATA			-6						// invalid data found in gl2 data file
#define GALREP_ECDATA_UNAVAILABLE		-7						// data for the requested curve and/or prime is not in the database
#define GALREP_GL2DATA_UNAVAILABLE		-8						// data for the group is not in the database
#define GALREP_SINGULAR_CURVE			-9						// the specified elliptic curve is singular
#define GALREP_OUT_OF_PRIMES			-10						// not enough primes in ec data to determine galois image within errbnd
#define GALREP_CCDATA_LOAD_FAILED		-11						// failed to load precomputed gl2 cc data
#define GALREP_CCSTATS_LOAD_FAILED		-12						// failed to load precomputed gl2 cc statistics
#define GALREP_BAD_CC_ID				-13						// bad cc id. must be in the range 0 to gl2_cc_count(ell) - 1
#define GALREP_BAD_CC_GEN_ID			-14						// bad cc generator id, must be in the range 0 to gl2_cc_gencnt(ell,id) - 1
#define GALREP_BAD_DET				-15						// determinant must be in the range 1 to ell-1, or -1 for wildcard
#define GALREP_BAD_TR					-16						// trace must be in the range 0 to ell-1 or -1 for wildcard
#define GALREP_MALLOC_FAILED			-17						// memory allocation error

#define GALREP_ECDATA_FILENAME			"galrep_ecdata.dat"			// default filename for precomputed ellipcitic curve data
#define GALREP_GL2DATA_FILENAME		"galrep_gl2data.dat"		// default filename for precomputed conjugacy class data for subgroups of GL(2,Z/mZ)
#define GALREP_DEFAULT_ERRBND			50						// error at most 2^{-100} (heuristically)
#define GALREP_MAX_ERRBND				1000					// error at most 2^{-1000} (heuristically)

#define GALREP_CC_FLAG_ABELIAN			1						// subgroup is abelian
#define GALREP_CC_FLAG_X_TRACE			2						// not every trace occurs
#define GALREP_CC_FLAG_X_N				4						// not every value of det+1-trace occurs
#define GALREP_CC_FLAG_E1				8						// 1 is an eignvalue of every element

int galrep_load (void);											// Loads default versions of precomputed data.  Called implicitly whenever needed (so no need to ever call this)
void galrep_unload (void);										// Unloads all precomputed data, freeing all dynamically allocated memory.
int galrep_ecdata_load (char *filename);							// Loads specified version of precomputed elliptic curve data (supercedes previously loaded data, if any)
int galrep_gl2data_load (char *filename);							// Loads specified version of precomputed gl2 conjugacy class data (supercedes previously loaded data, if any)

int galrep_ecdata_maxp (void);									// Returns largest prime p for which elliptic curve data has been loaded (0 if none)
int galrep_gl2data_maxl (void);									// Returns largest prime ell for which gl2 conjugacy class data has been loaded (0 if none)

// Elliptic curves over Q are specified in short Weierstrass fom y^2 = x^3 + Ax + B, with A and B integers (we never reduce mod 2 or 3)
// In the functions below errbnd must be zero (or equal to the default) this is reserved for future use
// A negative return value indicates an error, nonnegative indicates success

int galrep_ec_modl_image (int ell, mpz_t A, mpz_t B, int errbnd);					// returns cc id of the modl galois image (id 0 corresponds to the full group GL(2,Z/ellZ))
int galrep_ec_modl_images (int *ccs, int min, int max, mpz_t A, mpz_t B, int errbnd);	// sets ccs[i] for each ell in [min,max], with the first prime >= min in ccs[0], returns # of primes p used
int galrep_ec_non_surjective (int min, int max, mpz_t A, mpz_t B, int errbnd);		// returns a positive integer with the ith bit set iff the ith prime ell is in [min,max] and the modl image is nonsurjective

// long datatype versions of the mpz functions above
int galrep_ec_modl_image_i (int ell, long A, long B, int errbnd);
int galrep_ec_modl_images_i (int *ccs, int min, int max, long A, long B, int errbnd);
int galrep_ec_non_surjective_i (int min, int max, long A, long B, int errbnd);

int galrep_gl2_cc_count (int ell);									// number of subgroup conjugacy classes in database for GL(2,Z/ellZ) (only includes those that span all determinants)
int galrep_gl2_cc_index (int ell, int id);								// returns the index of the specified subgroup conjugacy class in GL(2,Z/ellZ)
int galrep_gl2_cc_order (int ell, int id);								// returns the order of the specified subgroup conjugacy class in GL(2,Z/ellZ)
int galrep_gl2_cc_tag (int ell, int id);								// obsolete, for internal use only
int galrep_gl2_cc_gen_count (int ell, int id);							// number of generators in a polycyclic+perfect presentation of a representative subgroup for the specified cc in GL(2,Z/ellZ)
int galrep_gl2_cc_gen (int ell, int id, int n);							// returns a positive integer encoding the nth generator A, defined as 2^24*A[1][1] + 2^16*A[1][0] + 2^8*A[0][1] + A[0][0]
int galrep_gl2_cc_freq (int ell, int id, int det, int tr);					// returns a count of the # of elements with specified det, tr.  Use -1 for wildcards.
int galrep_gl2_cc_flags (int ell, int id);								// returns a positive integer whose bits correspond to the flags defined by GALREP_CC_FLAG_* above

#endif
