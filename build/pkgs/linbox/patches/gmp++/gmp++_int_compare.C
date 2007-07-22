// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_compare.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_compare.C,v 1.5 2007-01-11 18:42:51 jgdumas Exp $
// ==========================================================================

#include "gmp++/gmp++.h"

// returns 1 if a > b, 0 if a == b and -1 otherwise.
int compare(const Integer &a, const Integer& b)
{
   return mpz_cmp ( (mpz_ptr)&a.gmp_rep, (mpz_ptr)&b.gmp_rep );
}

int absCompare(const Integer &a, const Integer &b)
{
   return mpz_cmpabs( (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep));
}

int Integer::operator != (const int l) const
{ return mpz_cmp_si ( (mpz_ptr)&gmp_rep, l ) != 0; }

int Integer::operator != (const long l) const
{ return mpz_cmp_si ( (mpz_ptr)&gmp_rep, l ) != 0; }

//unsigned long ops added by Dan Roche, 6-26-04
int Integer::operator != (const unsigned long l) const
{ return mpz_cmp_ui ( (mpz_ptr)&gmp_rep, l ) != 0; }

int Integer::operator > (const unsigned long l) const
{ return mpz_cmp_ui((mpz_ptr)&gmp_rep, l) > 0; }

int Integer::operator < (const unsigned long l) const
{ return mpz_cmp_ui((mpz_ptr)&gmp_rep, l) < 0; }

int Integer::operator > (const int l) const
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) > 0; }

int Integer::operator > (const long l) const
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) > 0; }

int Integer::operator < (const int l) const
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) < 0; }

int Integer::operator < (const long l) const
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) < 0; }

