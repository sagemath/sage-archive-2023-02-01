#ifndef __GMPplusplus_CSTOR_C__
#define __GMPplusplus_CSTOR_C__
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_cstor.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_cstor.C,v 1.3 2007-01-11 18:42:51 jgdumas Exp $
// ==========================================================================
#include <iostream>
#include "gmp++/gmp++.h"


//------------------------------------- predefined null and one
const Integer Integer::zero(0UL);
const Integer Integer::one(1UL);


// -- Integer(const char *s)
Integer::Integer(const char *s)
{
  mpz_init_set_str((mpz_ptr)&gmp_rep, s, 10);
}


Integer& Integer::copy(const Integer &n)
{
  if (this == &n) return *this;
  mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
  return *this ;
}

void importWords(Integer& x, size_t count, int order, int size, int endian, size_t nails, const void* op) {
  mpz_import( (mpz_ptr)&(x.gmp_rep), count, order, size, endian, nails, op);
}

#endif

