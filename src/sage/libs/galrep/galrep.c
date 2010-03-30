#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "gmp.h"
#include "galrep.h"


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

// NB: The int datatype is assumed to hold at least 32-bits

#define MASK_BITS			32
typedef uint32_t bitmask_t;

static int elltab[GALREP_ELL_COUNT] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59 };
static int ellitab[GALREP_MAX_ELL+1];		// ellitab[ell] = i iff elltab[i] = ell, set to -1 ow.  initialized by galrep_load_gl2data

/*
	General utility stuff for barebones modular arithmetic and bitmask manipulation.
	Note that we assume throughout that the modulus p is small enough so that 4*p^2 fits in an unsigned int
	We include this utility code here just to make things self contained.
*/
static inline void swab16 (uint16_t *x) { *x = ((*x&(uint16_t)0x00ffU) << 8) | ((*x&(uint16_t)0xff00U) >> 8); }
static inline void swab32 (uint32_t *x) { *x = ((*x&(uint32_t)0x000000ffU) << 24) | ((*x&(uint32_t)0x0000ff00U) << 8) | ((*x&(uint32_t)0x00ff0000U) >> 8) | ((*x&(uint32_t)0xff000000U) >> 24); }

double log2 (double x);

// the modular aritmetic below is about twice as fast as the % operator on an AMD Athlon 64 (YMMV)
static unsigned int ui_mod_i (unsigned int x, unsigned int m, double mi)
	{ register unsigned int t, z;  t = mi * x - 0.5;  z = x-t*m;  if ( z >= m ) z-= m; return z; }
static unsigned int ui_mod (unsigned int x, unsigned int m) { return ui_mod_i (x, m, 1.0/m); }
static int i_mod(int x, int m) { register int t;  if ( x >= 0 ) return ui_mod(x,m);  t = - (int)ui_mod((unsigned int)(-x),m);  return (t<0?t+m:t); }
int gcdext (int a, int b, int *x, int *y);
int legendre (int a, int b);
static int inverse_mod_p (int a, int p) { int y;  if ( gcdext(p, a, 0, &y) != 1 ) return 0; else return (y >= 0 ? y : y+p ); }
static int residue_mod_p (int a, int p) { return ( legendre (a, p) < 0 ? 0 : 1 ); }

static inline void bitset (bitmask_t *mask, int i) { mask[i/MASK_BITS] |= (1UL<<(i&(MASK_BITS-1))); }
static inline void bitclr (bitmask_t *mask, int i) { mask[i/MASK_BITS] &= ~(1UL<<(i&(MASK_BITS-1))); }
static inline int bittest (bitmask_t *mask, int i) { return ( (mask[i/MASK_BITS]>>(i&(MASK_BITS-1))) & 1UL ); }
static inline void maskand (bitmask_t *mask, bitmask_t *mask1, bitmask_t *mask2, int words) { register int i; for ( i = 0 ; i < words ; i++ ) mask[i] = mask1[i] & mask2[i]; }
static inline void maskor (bitmask_t *mask, bitmask_t *mask1, bitmask_t *mask2, int words) { register int i; for ( i = 0 ; i < words ; i++ ) mask[i] = mask1[i] | mask2[i]; }
static inline int maskclear (bitmask_t *mask, int words) { register int i; for ( i = 0 ; i < words ; i++ ) mask[i] = 0; }
static inline int masksetall (bitmask_t *mask, int words) { register int i; for ( i = 0 ; i < words ; i++ ) mask[i] = ~((bitmask_t)0UL); }
static inline int masknull (bitmask_t *mask, int words) { register int i; for ( i = 0 ; i < words ; i++ ) if ( mask[i] ) return 0; return 1; }
static inline int maskequal (bitmask_t *mask1, bitmask_t *mask2, int words) { register int i; for ( i = 0 ; i < words ; i++ )  if ( mask1[i] != mask2[i] ) return 0; return 1; }
static inline void maskcopy (bitmask_t *mask1, bitmask_t *mask2, int words) { register int i; for ( i = 0 ; i < words ; i++ ) mask1[i] = mask2[i]; }

/*
	ecdata routines go here - functions to load/unload/lookup precomputed group structure data
	for every elliptic curve over Fp for various small p > 3.
*/

/*
	Elliptic curve data for small primes in [4,32768] (the exact number depends on the file loaded).
	The jdata table includes for each j-invariant over F_p:
		1) the absolute value of the trace (stored in the high 11 bits)
		2) bit 0 indicating the sign of the trace for the twist of the form y^2=x^3+Ax+B for which A/B is a quadratic residue
		3) bits 1,2,3,4 indicating whether either twist has full 2,3,5,7-torsion mod p.
		    If bit 1 is set, both twists have 2-torsion, and for odd ell the twist for which p+1-t = 0 mod ell is the one with ell-torsion
	For ell from 11 to 59, if p is 1 mod ell then a p-terminated list of increasing j-invariants with full ell-torsion is stored.

	When p=1 mod 3, the data for j-invariant 0 should be ignored, since there are 6 twists, not 2, and we don't distinguish them (we could, but we don't)
	Similarly, when p=1 mod 4, the data for j-invariant 1728 should be ignored, since there are 4 twists, not 2.

	Note that we *could* distinguish the twists for 0 and 1728, we just don't *need* to, we will just ignore primes where this occurs.
	In the worst case, when we have a curve j-invariant 0 or 1728 over Q, we only use half the primes, but otherwise we use almost all of them.
*/

#define TORP			13
#define ECDATA_ID		12345

static struct ectab_entry_struct {
	uint16_t p;					// a prime p in the range [4,32768]
	uint16_t ptor;					// mask in which the nth bit is set iff p is 1 mod the (n+5)th prime (n=0 for 11, n=12 for 59)
	uint16_t *jdata;				// see details above
	uint16_t *torlist[TORP];			// torlist[i] points to a zero-terminated list of j-invariants with ell-torsion, where ell is the (n+5)th prime (n=0 for ell=11)
} *ecdata;
static int ecdatacnt;


void galrep_ecdata_unload (void)
{
	int i;

	if ( ecdata ) {
		for ( i = 0 ; i < ecdatacnt ; i++ ) if ( ecdata[i].jdata ) { free (ecdata[i].jdata); ecdata[i].jdata = 0; }
		free (ecdata);
	}
	ecdata = 0;  ecdatacnt = 0;
}

int galrep_ecdata_load (char *filename)
{
	FILE *fp;
	struct ectab_entry_struct *tabspace, *tabnext;
	char *s, *buf;
	uint16_t x16 , pcnt, *ps, *pe;
	int i, j, swab_flag, bytes;

	if ( ecdata ) galrep_ecdata_unload();
	if ( ! filename || ! filename[0] ) filename = (char*)GALREP_ECDATA_FILENAME;
	fp = fopen (filename, "rb");
	if ( ! fp ) return GALREP_FILE_NOT_FOUND;
	if ( fread(&x16,2,1,fp) != 1 ) goto read_error;
	if ( x16 != ECDATA_ID ) { swab_flag = 1;  swab16(&x16); if ( x16 != ECDATA_ID ) { fclose(fp); return GALREP_BAD_ECDATA; } } else { swab_flag = 0; }
	if ( fread(&x16,2,1,fp) != 1 ) goto read_error;  if ( swab_flag ) swab16(&x16);
	ecdatacnt = x16;  bytes = ecdatacnt * sizeof(*ecdata);
	ecdata = (struct ectab_entry_struct *) malloc(bytes);
	if ( ! ecdata ) { fclose(fp); galrep_ecdata_unload();  return GALREP_MALLOC_FAILED; }
	for ( i = 0 ; i < ecdatacnt ; i++ ) {
		if ( fread(&x16,2,1,fp) != 1 ) goto read_error;  if ( swab_flag ) swab16(&x16);  ecdata[i].p = x16;
		if ( fread(&x16,2,1,fp) != 1 ) goto read_error;  if ( swab_flag ) swab16(&x16);  ecdata[i].ptor = x16;
		if ( fread(&x16,2,1,fp) != 1 ) goto read_error;  if ( swab_flag ) swab16(&x16);  pcnt = x16;
		if ( pcnt < ecdata[i].p ) { fclose(fp); printf ("bad pcnt=%d for prime %d in ECDATA file %s\n", pcnt, ecdata[i].p, filename); galrep_ecdata_unload();  return GALREP_BAD_ECDATA; }
		ecdata[i].jdata = (uint16_t *) malloc(2*pcnt);
		pe = ecdata[i].jdata + pcnt;
		if ( ! ecdata[i].jdata ) { fclose (fp); galrep_ecdata_unload(); return GALREP_MALLOC_FAILED; }
		if ( fread(ecdata[i].jdata,2,pcnt,fp) != pcnt ) goto read_error;
		if ( swab_flag ) for ( ps = ecdata[i].jdata ; ps < pe ; ps++ ) swab16(ps);
		ps = ecdata[i].jdata + ecdata[i].p;
		for ( j = 0 ; j < TORP ; j++ ) {
			if ( ecdata[i].ptor & (((uint16_t)1)<<j) ) {
				ecdata[i].torlist[j] = ps;
				for ( ; ps < pe && *ps<ecdata[i].p ; ps++ );			// be careful not to read past the end of the buffer if it isn't properly terminated
				if ( ps++==pe ) { fclose(fp); printf ("invalid record format for prime %d in ECDATA file %s\n", ecdata[i].p, filename); galrep_ecdata_unload(); return GALREP_BAD_ECDATA; }
			} else {
				ecdata[i].torlist[j] = 0;
			}
		}
		if ( ps != pe ) { fclose(fp); printf ("invalid record format for prime %d in ECDATA file %s\n", ecdata[i].p, filename); galrep_ecdata_unload(); return GALREP_BAD_ECDATA; }
	}
	fclose(fp);
	/* printf("Loaded elliptic curve group data for %d primes up to %d from file %s\n", ecdatacnt, ecdata[ecdatacnt-1].p, filename);*/
	return 0;

read_error:
	fclose(fp);  galrep_ecdata_unload();
	return GALREP_FILE_READ_ERROR;
}

int galrep_ecdata_maxp (void) { return ( !ecdatacnt ? 0 : ecdata[ecdatacnt-1].p ); }

/*
	Given an elliptic curve E in the form y^2=x^3+Ax+B over F_p , returns the trace and optionally a bitmask
	identifying full ell-torsion subgroups (the ith bit is set iff for the (i+1)st prime ell either E or its twist has full ell-torsion over F_p).
	The prime p is specified by its index in ecdata, which is pi(p)-3 (i.e. pi=0 corresponds to 5).

	On input, the bitmask tormask should have the bits set for which torsion information is desired.
	Bits above 17 (corresponding to 59) will be ignored.
*/
int ecdata_lookup (int *trace, uint32_t *tormask, int A, int B, int pi)
{
	uint16_t *jp;
	uint32_t mask;
	register int i,j,k,t1,t2,jinv,p;									// we assume an int is at least 32 bits

	if ( ! ecdata || pi >= ecdatacnt ) return GALREP_ECDATA_UNAVAILABLE;
	p = ecdata[pi].p;

	if ( ! B ) {													// handle j-invariant 1728 separately
		if ( ! A ) return GALREP_SINGULAR_CURVE;					// this can happen when reducing mod p
		if ( !(p&2) ) return GALREP_ECDATA_UNAVAILABLE;			// ignore requests for j-invariant 1728 when p is 1 mod 4, we don't distinguish quartic twists
		*trace = 0;											// trace is always zero when p = 3 mod 4
		*tormask &= ( residue_mod_p (A,p) ? 0 : 1 );					// we can't have full odd-torsion, and we get full 2-torsion iff A is a non-residue (for p = 3 mod 4)
		return 0;
	}
	if ( ! A && (p%3)==1 ) return GALREP_ECDATA_UNAVAILABLE;		// ignore requests for j-invariant 0 when p is 1 mod 3, we don't distinguish sextic twists

	t1 = ui_mod(A*A,p); t1 = ui_mod(4*t1*A,p);						// t1 = 4A^3
	t2 = ui_mod(B*B,p);
	t2 = ui_mod(t1+27*t2,p);									// t2 = 4A^3 + 27B^2
	if ( ! t2 ) return GALREP_SINGULAR_CURVE;
	t2 = inverse_mod_p (t2, p);
	t1 = ui_mod(t1*t2, p);
	jinv = ui_mod(1728*t1,p);									// j = 1728 * 4A^3 / (4A^3+27B^2)
	t1 = ecdata[pi].jdata[jinv];									// lookup data for E/Fp in the database via its j-invariant
	*trace = t1 >> 5;											// get the absolute value of the trace
	if ( t1&0x10 ) *trace = - (*trace);								// check sign bit
	if ( A ) {
		if ( ! residue_mod_p (ui_mod(A*B,p),p)  ) *trace = -(*trace);		// adjust trace for the twist we have.  Note that B/A is a residue iff A*B is.
	} else {
		if ( ! residue_mod_p (B,p) ) *trace = -(*trace);				// when A is 0, we can just use B to distinguish twists (and there are only 2 because p=2 mod 3)
	}
	mask = (t1 & 0xf);											// get 2,3,5,7-torsion info
	*tormask &= (0xfffffff0|mask);									// set 2,3,5,7-torsion info
	mask = ecdata[pi].ptor;										// find out what other torsion info is possible for this p (i.e. ell for which p=1 mod ell)
	mask = *tormask & (mask << 4);								// restrict to ell's the caller actually cares about
	if ( ! mask ) { *tormask &= (0xf|mask); return 0;	}				// if nothing to do, return
	for ( i = 0 ; i < TORP ; i++ ) {									// for each ell from 11 to 59
		if ( mask&(((uint32_t)1)<<(i+4)) ) {						// if we are interested in this ell
			for ( jp = ecdata[pi].torlist[i] ; *jp < p && *jp != jinv ; jp++ );// do a linear search of the list of j-invariants with ell-torsion
			if ( *jp != jinv ) mask &= ~((uint32_t)1<<(i+4));			// if our j-invariant wasn't found, clear the corresponding bit
		}
	}
	*tormask &= (0xf|mask);
	return 0;
}

#define MAX_MASK_SIZE		(GALREP_MAX_CCS/32+((GALREP_MAX_CCS&0x1f)?1:0))	// largest possible mask size

struct cc_struct {
	int tag;									// legacy identifier, now deprecated
	int dups;									// number of distinct conjugacy classes with the same signature -- we only store data for one of these
	int order;								// subgroup size for this conjugacy class
	int pcnt;									// minimum number of primes we need to test to obtain a correct result with prob > 1-1/2^MAX_ERRBND (heuristically)
	uint32_t flags;								// bitmask identitying certain properties of the conjugacy class, see GALREP_CC_FLAG_xxx definitions in galrep.h
	int gencnt;								// number of generators for representative subgroup
	uint32_t gens[GALREP_MAX_GENERATORS];		// each generating matrix A is encoded as A[1][1]<<24 + A[1][0]<<16 + A[0][1]<<8 + A[0][0]
	bitmask_t upmask[MAX_MASK_SIZE];			// bit i is set for each class whose signature contains the signature of the current class
	uint16_t *counts;							// points to an array of counts -- the ith element counts non-trivial elements with i=(det-1)*ell+tr
};

static struct {
	int ell;									// currently the ith entry in this table always has ell equal to the (i+1)st prime (i.e. entry 0 has ell=2)
	int msize;								// size of subgroup masks for this ell, in words
	int cccnt;								// number of conjugacy classes in GL(2,Z/ellZ) that span determinants  with distinct signatures
	bitmask_t fullmask[MAX_MASK_SIZE];			// points to mask with a bit set for every proper subgroup of GL(2,Z/ellZ)
	bitmask_t *masks;							// array of ell*(ell-1) subgroup masks, ith entry corresponds to i=(det-1)*ell+tr
	struct cc_struct *ccs;						// points to and array of cc_count cc_structs
} gl2data[GALREP_ELL_COUNT];
static int gl2datacnt;							// number of entries currently loaded


#define GL2DATA_ID	123456789


void galrep_gl2data_unload (void)
{
	int i, j;

	for ( i = 0 ; i < gl2datacnt ; i++ ) {
		if ( gl2data[i].masks ) { free (gl2data[i].masks); gl2data[i].masks = 0; }
		if ( gl2data[i].ccs ) {
			for ( j = 0 ; j < gl2data[i].cccnt ; j++ ) if ( gl2data[i].ccs[j].counts ) free (gl2data[i].ccs[j].counts);
			gl2data[i].ccs = 0;
		}
	}
	gl2datacnt = 0;
}


int galrep_gl2data_load (char *filename)
{
	FILE *fp;
	uint8_t x8;
	uint16_t x16;
	uint32_t x32;
	int i, j, k, ell, msize, swab;

	// intialize inverse elltab
	for ( i = 0 ; i <= GALREP_MAX_ELL ; i++ ) ellitab[i] = -1;
	for ( i = 0 ; i < GALREP_ELL_COUNT ; i++ ) ellitab[elltab[i]] = i;
	if ( ! filename || ! filename[0] ) filename = (char*) GALREP_GL2DATA_FILENAME;
	fp = fopen (filename, "rb");
	if ( ! fp ) return GALREP_FILE_NOT_FOUND;
	if ( (i=fread (&x32,4,1,fp)) != 1 ) { printf ("fread returned %d\n", i); goto read_error; }
	if ( x32 != GL2DATA_ID ) { swab  = 1; swab32(&x32); if ( x32 != GL2DATA_ID ) { fclose(fp); return GALREP_BAD_GL2DATA; } } else swab = 0;

	for ( i = 0 ; i < GALREP_ELL_COUNT; i++ ) {
		if ( fread (&x8,1,1, fp) != 1 ) break;  gl2data[i].ell = ell = x8;
		if ( gl2data[i].ell != elltab[i] ) { fclose(fp); return GALREP_BAD_GL2DATA; }
		if ( fread (&x8,1,1, fp) != 1 ) goto read_error;  gl2data[i].msize = msize = x8;
		if ( fread (&x16,2,1, fp) != 1 ) goto read_error;  if ( swab ) swab16(&x16);  gl2data[i].cccnt = x16;
		if ( gl2data[i].cccnt > GALREP_MAX_CCS || 32*gl2data[i].msize < gl2data[i].cccnt ) { fclose(fp); galrep_gl2data_unload (); return GALREP_BAD_GL2DATA; }
		if ( fread (gl2data[i].fullmask,4,msize,fp) != msize ) goto read_error;
		if ( swab ) for ( k = 0 ; k < msize ; k++ ) swab32(gl2data[i].fullmask+k);
		gl2data[i].masks = (bitmask_t *) malloc (ell*(ell-1)*msize*sizeof(bitmask_t));
		if ( ! gl2data[i].masks ) { fclose(fp); galrep_gl2data_unload (); return GALREP_MALLOC_FAILED; }
		if ( fread(gl2data[i].masks,4,ell*(ell-1)*msize,fp) != ell*(ell-1)*msize ) goto read_error;
		if ( swab ) for ( j = 0 ; j < ell*(ell-1)*msize ; j++ ) swab32(gl2data[i].masks+j);
		gl2data[i].ccs = (struct cc_struct *) calloc(gl2data[i].cccnt,sizeof(struct cc_struct));
		if ( ! gl2data[i].ccs ) { fclose(fp); galrep_gl2data_unload (); return GALREP_MALLOC_FAILED; }
		for ( j = 0 ; j < gl2data[i].cccnt ; j++ ) {
			if ( fread (&x16,2,1,fp) != 1 ) goto read_error;  if ( swab ) swab16(&x16);  	gl2data[i].ccs[j].tag = x16;
			if ( fread (&x8,1,1,fp) != 1 ) goto read_error;  gl2data[i].ccs[j].dups = x8;
			if ( fread (&x8,1,1,fp) != 1 ) goto read_error;  gl2data[i].ccs[j].gencnt = x8;
			if ( gl2data[i].ccs[j].gencnt > GALREP_MAX_GENERATORS ) { fclose(fp); return GALREP_BAD_GL2DATA; }
			if ( fread (&x32,4,1,fp) != 1 ) goto read_error;  if ( swab ) swab32(&x32);  	gl2data[i].ccs[j].order = x32;
			if ( fread (&x32,4,1,fp) != 1 ) goto read_error;  if ( swab ) swab32(&x32);  	gl2data[i].ccs[j].flags = x32;
			if ( fread (&x32,4,1,fp) != 1 ) goto read_error;  if ( swab ) swab32(&x32);  	gl2data[i].ccs[j].pcnt = x32;
			if ( fread (gl2data[i].ccs[j].upmask,4,msize,fp) != msize ) goto read_error;
			if ( swab ) for ( k = 0 ; k < msize ; k++ ) swab32(gl2data[i].ccs[j].upmask+k);
			if ( fread (gl2data[i].ccs[j].gens,4,gl2data[i].ccs[j].gencnt,fp) != gl2data[i].ccs[j].gencnt ) goto read_error;
			if ( swab ) for ( k = 0 ; k < gl2data[i].ccs[j].gencnt ; k++ ) swab32(gl2data[i].ccs[j].gens+k);
			gl2data[i].ccs[j].counts = (uint16_t *) malloc(ell*(ell-1)*sizeof(uint16_t));
			if ( ! gl2data[i].ccs[j].counts ) { fclose (fp); galrep_gl2data_unload (); return GALREP_MALLOC_FAILED; }
			if ( fread (gl2data[i].ccs[j].counts,2,ell*(ell-1),fp) != ell*(ell-1) ) goto read_error;
			if ( swab ) for ( k = 0 ; k < ell*(ell-1) ; k++ ) swab16(gl2data[i].ccs[j].counts+k);
		}
	}
	fclose (fp);
	gl2datacnt = i;
	/* printf("Loaded conjugacy class data for GL(2,Z/ellZ) for %d primes up to %d from file %s\n", gl2datacnt, gl2data[i-1].ell, filename);*/
	return 0;
read_error:
	fclose(fp);  galrep_ecdata_unload();
	return GALREP_FILE_READ_ERROR;
}

int galrep_gl2data_maxl (void) { return ( !gl2datacnt ? 0 : gl2data[gl2datacnt-1].ell ); }

int galrep_gl2_cc_count (int ell)
{
	int i;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	return gl2data[i].cccnt;
}

int galrep_gl2_cc_order (int ell, int id)
{
	int i;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	if ( id < 0 || id >= gl2data[i].cccnt ) return GALREP_BAD_CC_ID;
	return gl2data[i].ccs[id].order;
}

int galrep_gl2_cc_index (int ell, int id)
{
	int n, o;

	o = galrep_gl2_cc_order(ell,id);
	if ( o < 0 ) return o;
	n = ell*(ell-1)*(ell-1)*(ell+1);		// we assume ell is small enough to avoid overflow here
	return n/o;
}

int galrep_gl2_cc_tag (int ell, int id)
{
	int i;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	if ( id < 0 || id >= gl2data[i].cccnt ) return GALREP_BAD_CC_ID;
	return gl2data[i].ccs[id].tag;
}

int galrep_gl2_cc_dups (int ell, int id)
{
	int i;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	if ( id < 0 || id >= gl2data[i].cccnt ) return GALREP_BAD_CC_ID;
	return gl2data[i].ccs[id].dups;
}

int galrep_gl2_cc_flags (int ell, int id)
{
	int i;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	if ( id < 0 || id >= gl2data[i].cccnt ) return GALREP_BAD_CC_ID;
	return gl2data[i].ccs[id].flags;
}

int galrep_gl2_cc_gen_count (int ell, int id)
{
	int i;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	if ( id < 0 || id >= gl2data[i].cccnt ) return GALREP_BAD_CC_ID;
	return gl2data[i].ccs[id].gencnt;
}

int galrep_gl2_cc_gen (int ell, int id, int n)
{
	int i;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	if ( id < 0 || id >= gl2data[i].cccnt ) return GALREP_BAD_CC_ID;
	if ( n < 0 || n > gl2data[i].ccs[id].gencnt ) return GALREP_BAD_CC_GEN_ID;
	return gl2data[i].ccs[id].gens[n];
}

int galrep_gl2_cc_freq (int ell, int id, int det, int tr)
{
	int i,k,cnt;

	if ( ell < 0 || ell > GALREP_MAX_ELL || (i=ellitab[ell]) < 0 ) return GALREP_BAD_ELL;
	if ( i >= gl2datacnt ) return GALREP_GL2DATA_UNAVAILABLE;
	if ( id < 0 || id >= gl2data[i].cccnt ) return GALREP_BAD_CC_ID;
	if ( det < -1 || det == 0 || det >= ell ) return GALREP_BAD_DET;
	if ( tr < -1 || tr >= ell ) return GALREP_BAD_TR;
	if ( det >= 1 && tr >= 0 ) {
		k = (det-1)*ell+tr;
		cnt = gl2data[i].ccs[id].counts[k];
		if ( det==1 && (tr==2 || (ell==2&&!tr)) ) cnt++;
	} else if ( tr >= 0 ) {
		if ( tr==2 || (!tr && ell==2) ) cnt = 1; else cnt = 0;
		for ( det = 1 ; det < ell ; det++ ) {
			k = (det-1)*ell+tr;
			cnt += gl2data[i].ccs[id].counts[k];
		}
	} else if ( det >= 0 ) {
		if ( det==1 ) cnt = 1; else cnt = 0;
		for ( tr = 0 ; tr < ell ; tr++ ) {
			k = (det-1)*ell+tr;
			cnt += gl2data[i].ccs[id].counts[k];
		}
	} else {
		cnt = gl2data[i].ccs[id].order;
	}
	return cnt;
}

struct curve_ctx {
	int start_index, end_index;
	bitmask_t submask[GALREP_ELL_COUNT][MAX_MASK_SIZE];
	int pcnt[GALREP_ELL_COUNT];
	int done[GALREP_ELL_COUNT];
	int cc[GALREP_ELL_COUNT];
	int done_count;
	int errbnd;
};

void process_frobenius (struct curve_ctx *ctx, int p, int ap, uint32_t tormask);
int setup_ctx (struct curve_ctx *ctx, int minell, int maxell, int errbnd);

int galrep_ec_modl_image (int ell, mpz_t A, mpz_t B, int errbnd)
	{ int cc, err;  err = galrep_ec_modl_images (&cc, ell, ell, A, B, errbnd);  if ( err < 0 ) return err;  return cc; }

int galrep_ec_modl_images (int *ccs, int min, int max, mpz_t A, mpz_t B, int errbnd)
{
	uint32_t tormask;
	int a, b, p, ap;
	struct curve_ctx ctx;
	int ell_count, err;
	register int i,j;

	if ( (err=setup_ctx (&ctx, min, max, errbnd)) < 0 ) return err;
	ell_count = ctx.end_index - ctx.start_index;
	for ( i = 0 ; i < ecdatacnt ; i++ ) {
		p = ecdata[i].p;
		a = (int) mpz_fdiv_ui (A, p);  b = (int) mpz_fdiv_ui (B, p);
		tormask = 0x1FF;
		if ( ecdata_lookup (&ap, &tormask, a, b, i) < 0 ) continue;		// skip primes where we get an error, e.g. primes of bad reduction
		process_frobenius (&ctx, p, ap, tormask);
		if ( ctx.done_count == ell_count ) break;
	}
	if ( i==ecdatacnt ) return GALREP_OUT_OF_PRIMES;
	for ( j = 0 ; j < ell_count ; j++ ) ccs[j] = ctx.cc[ctx.start_index+j];
	return i+1;											// return the number of primes we used, just in case someone is curious
}

int galrep_ec_non_surjective (int min, int max, mpz_t A, mpz_t B, int errbnd)
{
	int ccs[GALREP_ELL_COUNT];
	int i, j, err, retval;

	if ( max > GALREP_MAX_ELL ) return GALREP_BAD_ELL;
	err = galrep_ec_modl_images (ccs, min, max, A, B, errbnd);
	if ( err < 0 ) return err;
	retval = 0;
	for ( i = 0 ; elltab[i] < min ; i++ );
	for ( j = i ; elltab[i] <= max ; i++) if ( ccs[i-j] > 0 ) retval |= (1<<i);
	return retval;
}

int galrep_ec_modl_image_i (int ell, long A, long B, int errbnd)
	{ int cc, err;  err = galrep_ec_modl_images_i (&cc, ell, ell, A, B, errbnd);  if ( err < 0 ) return err;  return cc; }

int galrep_ec_modl_images_i (int *ccs, int min, int max, long A, long B, int errbnd)
{
	uint32_t tormask;
	int a, b, p, ap;
	struct curve_ctx ctx;
	int ell_count, err;
	register int i,j;

	if ( (err=setup_ctx (&ctx, min, max, errbnd)) < 0 ) return err;
	ell_count = ctx.end_index - ctx.start_index;
	for ( i = 0 ; i < ecdatacnt ; i++ ) {
		p = ecdata[i].p;
		a = i_mod (A, p);  b = i_mod (B, p);
		tormask = 0x1FF;
		if ( ecdata_lookup (&ap, &tormask, a, b, i) < 0 ) continue;		// skip primes where we get an error, e.g. primes of bad reduction
		process_frobenius (&ctx, p, ap, tormask);
		if ( ctx.done_count == ell_count ) break;
	}
	if ( i==ecdatacnt ) return GALREP_OUT_OF_PRIMES;
	for ( j = 0 ; j < ell_count ; j++ ) ccs[j] = ctx.cc[ctx.start_index+j];
	return i+1;											// return the number of primes we used, just in case someone is curious
}

int galrep_ec_non_surjective_i (int min, int max, long A, long B, int errbnd)
{
	int ccs[GALREP_ELL_COUNT];
	int i, j, err, retval;

	if ( max > GALREP_MAX_ELL ) return GALREP_BAD_ELL;
	err = galrep_ec_modl_images_i (ccs, min, max, A, B, errbnd);
	if ( err < 0 ) return err;
	retval = 0;
	for ( i = 0 ; elltab[i] < min ; i++ );
	for ( j = i ; elltab[i] <= max ; i++) if ( ccs[i-j] > 0 ) retval |= (1<<i);
	return retval;
}

int setup_ctx (struct curve_ctx *ctx, int minell, int maxell, int errbnd)
{
	int i, err;

	if ( maxell > GALREP_MAX_ELL ) return GALREP_BAD_ELL;
	if ( errbnd < 0 || errbnd > GALREP_MAX_ERRBND ) return GALREP_BAD_ERRBND;
	if ( ! errbnd ) errbnd = GALREP_DEFAULT_ERRBND;
	if ( ! ecdatacnt && (err=galrep_ecdata_load(0)) ) return err;
	if ( ! gl2datacnt && (err=galrep_gl2data_load(0)) ) return err;

	ctx->done_count = 0;  ctx->errbnd = errbnd;
	for ( i = 0 ; i < GALREP_ELL_COUNT && elltab[i] < minell ; i++ );
	ctx->start_index = i;
	while ( i < GALREP_ELL_COUNT && elltab[i] <= maxell ) i++;
	ctx->end_index = i;
	if ( ctx->start_index == ctx->end_index )  return GALREP_BAD_ELL;
	for ( i = ctx->start_index ; i < ctx->end_index ; i++ ) { ctx->done[i] = ctx->pcnt[i] = 0;  maskcopy (ctx->submask[i], gl2data[i].fullmask, gl2data[i].msize); }
	return 0;
}

void process_frobenius (struct curve_ctx *ctx, int p, int ap, uint32_t tormask)
{
	register int i, j, k, ell, msize;

//	printf ("process_frobenius: p=%d, ap=%d, tor=%x (hex)\n", p, ap, tormask);
	for ( i = ctx->start_index ; i < ctx->end_index ; i++ ) {
		if  ( ctx->done[i] ) continue;
		ell = gl2data[i].ell;
		if ( p==ell ) continue;
		msize = gl2data[i].msize;
		ctx->pcnt[i]++;
		// Only process frobenius elements pi_p that act non-trivially mod ell, i.e. for which we do not have full ell-torsion mod p
		// Note that tormask doesn't distinguish the trace, so if ell is odd we also check that the group order is divisible by ell
		if ( ! (tormask&((uint32_t)1<<i)) || (i && ui_mod(p+1-ap,ell)) ) {
			j = ui_mod(p,ell);  k = i_mod(ap,ell);
			k = (j-1)*ell+k;
			maskand(ctx->submask[i],ctx->submask[i],gl2data[i].masks+k*msize,msize);
			if ( masknull(ctx->submask[i], msize) ) { ctx->cc[i] = 0; ctx->done[i] = 1; ctx->done_count++; continue; }
		}
		if ( ctx->pcnt[i] > ctx->errbnd && !(ctx->pcnt[i]&0x7) ) {	// only check every 8th prime to avoid a lot of fruitless checking
			// we do a linear search here, we could use a binary search or a hash table, but the lists are pretty short (< 200)
			for ( j = 0 ; j < gl2data[i].cccnt ; j++ ) if ( maskequal (ctx->submask[i],gl2data[i].ccs[j].upmask,msize) ) break;
			if ( j < gl2data[i].cccnt && ctx->pcnt[i]*GALREP_MAX_ERRBND >= gl2data[i].ccs[j].pcnt*ctx->errbnd ) { ctx->cc[i] = j; ctx->done[i] = 1; ctx->done_count++; continue; }
		}
	}
}


// both a and b must be nonnegative
int gcdext (int a, int b, int *x, int *y)
{
	register int q, r, s, t, r0, r1, s0, s1, t0, t1;

	if ( a < b ) return gcdext (b, a, y, x);
	if ( b == 0 ) {
		if ( x ) *x = 1;
		if ( y ) *y = 0;
		return a;
	}
	if ( x ) { s0 = 1;  s1 = 0; }
	if ( y ) { t1 = 1;  t0 = 0; }
	r0 = a;  r1 = b;
	while ( r1 > 0 ) {
		q = r0/r1;  r = r0 - q*r1;
		r0 = r1;  r1 = r;
		if ( y ) { t = t0 - q*t1;  t0 = t1;  t1 = t; }
		if ( x ) { s = s0 - q*s1;  s0 = s1;  s1 = s; }
	}
	if ( x ) *x = s0;
	if ( y ) *y = t0;
	return r0;
}

// both a and b must be nonnegative and b must be odd.  this is not verified
int legendre (int a, int b)
{
	register int r, k, v;

	k = 1;
	while ( a ) {
		for ( v = 0 ; ! (a&0x1) ; v++ ) a >>= 1;
		if ( v&0x1 )  if ( (b&0x7) ==3 || (b&0x7) ==5 ) k = -k;
		if ( (a&0x3) == 3 && (b&0x3) == 3 ) k = -k;
		r = a;   a = ui_mod(b,r);  b = r;
	}
	return ( b == 1 ? k : 0 );
}
