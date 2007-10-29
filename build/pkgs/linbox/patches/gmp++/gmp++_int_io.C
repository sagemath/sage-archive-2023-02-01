// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_io.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_io.C,v 1.5 2007-01-11 18:42:51 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <stdlib.h>
#include "gmp++/gmp++.h"

// Sortie nonsignee : 321321 meme si n = -321321, par exemple
std::ostream& absOutput(std::ostream &o, const Integer&n)
{
  int base = 10;

  unsigned long strSize = mpz_sizeinbase((mpz_ptr)&(n.gmp_rep), base) + 2;
  char *str = ::new char[strSize];
  mpz_get_str(str, base, (mpz_ptr)&(n.gmp_rep));
  if (sign(n) < 0) {
      char *str1 = &(str[1]) ;
      o << str1;
  }
  else o << str;
  delete [] str ;
  return o;
}

// Sortie signee : +321321 ou -321321, par exemple
std::ostream& Integer::print(std::ostream &o) const
{
  int base = 10;
  unsigned long strSize = mpz_sizeinbase((mpz_ptr)&(gmp_rep), base) + 2;
  char *str = new char[strSize];
  mpz_get_str(str, base, (mpz_ptr)&(gmp_rep));
// JGD 08.11.1999 : temporaire
//   if (sign(*this) > 0) o << '+' ;
  o << str;
  delete [] str ;
  return o;
}

Integer::operator std::string () const {
    std::string s;
    unsigned long strSize = mpz_sizeinbase((mpz_ptr)&(gmp_rep), 10) + 2;
    char *str = new char[strSize + 2];
    mpz_get_str(str, 10, (mpz_ptr)&(gmp_rep));
    s = std::string(str);
    delete [] str ;
    return s;
}



Integer::Integer(const std::vector<mp_limb_t>& v) {
 	size_t s = v.size();
	if (s) {
	 	mpz_init_set_ui((mpz_ptr)&gmp_rep, v[0]);
		Integer base(256), prod, tmp;
		prod = base = pow(base, (unsigned long)sizeof(mp_limb_t) );

		std::vector<mp_limb_t>::const_iterator vi = v.begin();
		for(++vi;vi != v.end();++vi) {
			mpz_mul_ui( (mpz_ptr)&tmp.gmp_rep, (mpz_ptr)&prod.gmp_rep, *vi);
	 		*this += tmp;
			prod *= base;
		}
	} else
	 	mpz_init( (mpz_ptr)&gmp_rep );

}

Integer::operator std::vector<mp_limb_t> () const {
        size_t s = mpz_size( (mpz_ptr)&(gmp_rep) );
        std::vector<mp_limb_t> v(s);
        std::vector<mp_limb_t>::iterator vi = v.begin();
        for(mp_size_t i = 0;vi != v.end();++vi, ++i) *vi = mpz_getlimbn( (mpz_ptr)& (gmp_rep) ,i);
        return v;
}



  // Entree au format de la sortie
std::istream& operator>> (std::istream& in, Integer& a)
{
   static long base[10] = {
     10,
     100,
     1000,
     10000,
     100000,
     1000000,
     10000000,
     100000000,
     1000000000
   } ;
   if (!in) return in ;
   // eat white
   in >> std::ws  ;

   // Base : 10^9, we read by packet of length 9
   // the char.
   char Tmp[10] ;
   int counter = 0 ;

   // Set the returned integer
   a = 0L ;
   char ch ;
   int sign = 1 ;

   // find a sign:
   in.get(ch) ;
   if ((ch != '+') && (ch != '-') && !((ch >= '0') && (ch <= '9')))
   {
      std::cerr << "Bad integer format: found: "<< ch ;
      std::cerr << ", in place of '+' '-' or a digit"<< std::endl ;
      return in ;
   }
   switch (ch) {
     case '+' : break ;
     case '-' : sign = -1 ; break ;
     default  : in.putback(ch) ; break ;
   }
   // eat white
   in >> std::ws  ;

   int noend = 1 ;
   while (noend)
   {
      counter = 0 ;

      // Read 9 digits or less
      while ((noend) && (counter < 9)) {
         in.get(ch) ;
         if (in.eof()) { noend = 0 ; }
         else if ((ch >= '0') && (ch <= '9')) Tmp[counter++] = ch ;
         else { noend = 0 ;  in.putback(ch) ; }
      }
      if (counter >0) {
         long l ;
         Tmp[counter] = '\0' ; // terminate the string
         l = atol(Tmp) ;
         a = a * base[counter-1] + l ;
      }
   }
   if (sign == -1) a = -a ;
   return in ;
}
