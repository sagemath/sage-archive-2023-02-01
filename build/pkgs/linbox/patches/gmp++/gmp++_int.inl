// ========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int.inl,v 1.8 2007-01-11 18:42:51 jgdumas Exp $
// ========================================================================
// Description:

#define GMP__ABS(l)     ((l) <0 ? -l : l)
#define GMP__SGN(l)    ((l) <0 ? -1 : (l >0 ? 1 : 0))


//-----------------------------~Integer()
inline Integer::~Integer() {  mpz_clear((mpz_ptr)&gmp_rep) ; }

//-------------------------------Integer(const Integer &n)
inline Integer::Integer(const Integer &n) {
    mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
}

//------------------------------------------operator = (const Integer &n)
inline Integer& Integer::logcpy(const Integer &n)
{
  if (this == &n) return *this;
  mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
  return *this;
}

// same as logcopy
inline Integer& Integer::operator = (const Integer &n) { return logcpy(n) ; }

//-----------------------------Integer(int n)
inline Integer::Integer(int n) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(uint n)
inline Integer::Integer(unsigned char n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(uint n)
inline Integer::Integer(unsigned int n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(long n)
inline Integer::Integer(long n) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(unsigned long n)
inline Integer::Integer(unsigned long n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
#include <stdio.h>
//-----------------------------Integer(long long n)
// log[10](2^8) < 2.408239966
inline Integer::Integer(long long n) {
char * tmp = new char[long(2.408239966*sizeof(long long))+1]; sprintf(tmp,"%lld",n);
 mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
delete [] tmp;
}

//-----------------------------Integer(unsigned long long n)
// log[10](2^8) < 2.408239966
inline Integer::Integer(unsigned long long n) {
char * tmp = new char[ long(2.408239966*sizeof(unsigned long long))+1];
sprintf(tmp,"%llu",n);
mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
delete [] tmp;
}
#endif


//-----------------------------Integer(double)
inline Integer::Integer(double d) { mpz_init_set_d((mpz_ptr)&gmp_rep, d) ; }


//-----------------------------Integer(const neutral n), default n = zero
/* Neutral is causing a problem
inline Integer::Integer(const Neutral n) {
  if (n == Neutral::zero) mpz_init_set_ui((mpz_ptr)&gmp_rep, 0L) ;
  else  mpz_init_set_ui((mpz_ptr)&gmp_rep, 1L) ;
}
*/

//-------------------------------------------------inline comparaison operators
inline int operator != (const Integer& a , const Integer& b)
  { return compare(a,b) != 0; }

inline int operator != (int l, const Integer& n)
  { return n.operator != (l); }

inline int operator != (long l, const Integer& n)
  { return n.operator != (l); }

inline int operator != (unsigned long l, const Integer& n)
  { return n.operator != (l); }

inline int operator == (const Integer& a, const Integer& b)
  {  return compare(a,b) == 0; }

inline int operator == (int l, const Integer& n)
  { return (! (n.operator != (l))); }

inline int operator == (long l, const Integer& n)
  { return (! (n.operator != (l))); }

inline int operator == (unsigned long l, const Integer& n)
  { return (! (n.operator != (l))); }

inline int operator == (const Integer& n, unsigned long l)
  { return (! (n.operator != (l))); }

inline int operator == (const Integer& n, int l)
  { return (! (n.operator != (l))); }

inline int operator == (const Integer& n, long l)
  { return (! (n.operator != (l))); }

inline int operator < (const Integer& a , const Integer& b)
  { return compare(a,b) < 0; }

inline int operator < (const int l, const Integer& n)
  { return n > l; }

inline int operator < (const long l, const Integer& n)
  { return n > l; }

inline int operator < (const unsigned long l, const Integer& n)
  { return n > l; }

inline int operator <= (const Integer& n, unsigned long l)
  {  return (! (n > l) ); }

inline int operator <= (unsigned long l, const Integer& n)
  {  return (! (n < l) );}

inline int operator >= (unsigned long l, const Integer& n)
  {  return (! (n < l) );}

inline int operator >= (const Integer& n, unsigned long l)
  {  return (! (n < l) );}

inline int operator > (int l, const Integer& n)
  { return n < l; }

inline int operator > (long l, const Integer& n)
  { return n < l; }

inline int operator > (unsigned long l, const Integer& n)
  { return n < l; }

inline int operator >  (const Integer& a , const Integer& b)
  { return compare(a,b) > 0; }

inline int operator <= (const Integer& a, const Integer& b)
  { return compare(a,b) <= 0; }

inline int operator <= (const Integer& n, int l)
  {  return (! (n > l) ); }

inline int operator <= (const Integer& n, long l)
  {  return (! (n > l) ); }

inline int operator <= (int l, const Integer& n)
  {  return (! (n < l) );}

inline int operator <= (long l, const Integer& n)
  {  return (! (n < l) );}

inline int operator >= (const Integer& a, const Integer& b)
  { return compare(a,b) >= 0; }

inline int operator >= (int l, const Integer& n)
  {  return (! (n > l) );}

inline int operator >= (long l, const Integer& n)
  {  return (! (n > l) );}

inline int operator >= (const Integer& n, int l)
  {  return (! (n < l) );}

inline int operator >= (const Integer& n, long l)
  {  return (! (n < l) );}


//----------------------------------arithmetic inline operators
inline Integer Integer::operator - () const
{
// JGD 18.06.1999
    Integer Res ;
    mpz_neg((mpz_ptr)&Res.gmp_rep, (mpz_ptr)&gmp_rep );
    return Res ;
}

// -- operator +
inline Integer operator + (const int l, const Integer& n) { return n + (long)l; }
inline Integer operator + (const unsigned int l, const Integer& n) { return n + (unsigned long)l; }
inline Integer operator + (const long l, const Integer& n) { return n + l; }
inline Integer operator + (const unsigned long l, const Integer& n) { return n + l; }
inline Integer operator + (const Integer& n, const int l) { return n + (long)l; }
inline Integer operator + (const Integer& n, const unsigned int l) { return n + (unsigned long)l; }

inline Integer& operator += (Integer& n, const int l) { return n += (long)l; }
inline Integer& operator += (Integer& n, const unsigned int l) { return n += (unsigned long)l; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
  inline Integer operator + (const Integer& n, const long long l) {return n + (Integer)l; }
  inline Integer operator + (const Integer& n, const unsigned long long l) {return n + (Integer)l; }
  inline Integer operator + (const long long l, const Integer& n) {return n+l;}
  inline Integer operator + (const unsigned long long l, const Integer& n) {return n+l;}
  inline Integer& operator += (Integer& n, const long long l) { return n += (Integer)l; }
  inline Integer& operator += (Integer& n, const unsigned long long l) { return n += (Integer)l; }
#endif


// -- operator -
inline Integer operator - (const int l, const Integer& n) { return -(n - (long)l); }
inline Integer operator - (const unsigned int l, const Integer& n) { return -(n - (unsigned long)l); }
inline Integer operator - (const long l, const Integer& n) { return -(n - l); }
inline Integer operator - (const unsigned long l, const Integer& n) { return -(n - l); }
inline Integer operator - (const Integer& n, const int l) { return n - (long)l; }
inline Integer operator - (const Integer& n, const unsigned int l) { return n - (unsigned long)l; }

inline Integer& operator -= (Integer& n, const int l) { return n -= (long)l; }
inline Integer& operator -= (Integer& n, const unsigned int l) { return n -= (unsigned long)l; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
  inline Integer operator - (const Integer& n, const long long l) {return n - (Integer)l; }
  inline Integer operator - (const Integer& n, const unsigned long long l) {return n - (Integer)l; }
  inline Integer operator - (const long long l, const Integer& n) {return n-l;}
  inline Integer operator - (const unsigned long long l, const Integer& n) {return n-l;}
  inline Integer& operator -= (Integer& n, const long long l) { return n -= (Integer)l; }
  inline Integer& operator -= (Integer& n, const unsigned long long l) { return n -= (Integer)l; }
#endif

// -- operator *
inline Integer operator * (const int l, const Integer& n) { return n * (long)l; }
inline Integer operator * (const unsigned int l, const Integer& n) { return n * (unsigned long)l; }
inline Integer operator * (const long l, const Integer& n) { return n * l; }
inline Integer operator * (const unsigned long l, const Integer& n) { return n * l; }
inline Integer operator * (const Integer& n, const int l) { return n * (long)l; }
inline Integer operator * (const Integer& n, const unsigned int l) { return n * (unsigned long)l; }

inline Integer& operator *= (Integer& n, const int l) { return n *= (long)l; }
inline Integer& operator *= (Integer& n, const unsigned int l) { return n *= (unsigned long)l; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
  inline Integer operator * (const Integer& n, const long long l) {return n * (Integer)l; }
  inline Integer operator * (const Integer& n, const unsigned long long l) {return n * (Integer)l; }
  inline Integer operator * (const long long l, const Integer& n) {return n*l;}
  inline Integer operator * (const unsigned long long l, const Integer& n) {return n*l;}
  inline Integer& operator *= (Integer& n, const long long l) { return n *= (Integer)l; }
  inline Integer& operator *= (Integer& n, const unsigned long long l) { return n *= (Integer)l; }
#endif

// -- operator /
inline Integer operator / (const int l, const Integer& n) { return Integer(l)/n; }
inline Integer operator / (const long l, const Integer& n) { return Integer(l)/n; }
inline Integer operator / (const Integer& n, const int l) { return n / (long)l; }
inline Integer operator / (const Integer& n, const unsigned int l) { return n / (unsigned long)l; }

inline Integer& operator /= (Integer& n, const int l) { if (l>=0) return n /= (unsigned long)l; else return  n = -(n / (unsigned long)-l); }
inline Integer& operator /= (Integer& n, const long l) { return n /= (unsigned long)l; }
inline Integer& operator /= (Integer& n, const unsigned int l) { return n /= (unsigned long)l; }

// -- operator %
inline Integer operator % (const int l, const Integer& n) { return Integer(l) % n; }
inline Integer operator % (const long l, const Integer& n) { return Integer(l) % n; }
inline Integer operator % (const Integer& n, const int l) { return n % (long)l; }
inline Integer operator % (const Integer& n, const unsigned int l) { return n % (unsigned long)l; }

inline Integer& operator %= (Integer& n, const int l) { return n %= (long)l; }
inline Integer& operator %= (Integer& n, const unsigned int l) { return n %= (unsigned long)l; }


//----------miscellaneous inline functions

inline int Integer::priv_sign() const { return mpz_sgn( (mpz_ptr)&gmp_rep ); }

inline int isOne(const Integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 1UL); }

inline int isZero(const Integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 0UL); }

inline int isZero(const short int a) { return a ==0; }
inline int isZero(const int a) { return a ==0; }
inline int isZero(const long a) { return a ==0; }
inline int isZero(const unsigned short int a) { return a ==0; }
inline int isZero(const unsigned int a) { return a ==0; }
inline int isZero(const unsigned long a) { return a ==0UL; }
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
inline int isZero(const unsigned long long a) { return a ==0ULL; }
inline int isZero(const long long a) { return a ==0LL; }
#endif

inline int sign(const Integer& a) { return a.priv_sign(); }

inline unsigned long length(const Integer& a) { return mpz_size( (mpz_ptr)&(a.gmp_rep) ) * sizeof(unsigned long); }

inline Integer abs(const Integer &n) { if (sign(n) >= 0) return n; return -n; }

inline size_t Integer::size() const { return  mpz_size( (mpz_ptr)&gmp_rep ) ; }

inline size_t Integer::size_in_base(int BASE) const { return  mpz_sizeinbase ((mpz_ptr)&gmp_rep, BASE);}

inline size_t Integer::bitsize() const { return  mpz_sizeinbase ((mpz_ptr)&gmp_rep, 2);}

inline unsigned long Integer::operator[](size_t i) const
{ if ( mpz_size( (mpz_ptr)&gmp_rep ) > i)
    return mpz_getlimbn( (mpz_ptr)&gmp_rep, i);
 else
     return 0;
}

//-------------------------------------------------inline >> & << operators
inline std::ostream& operator<< (std::ostream& o, const Integer& a) { return a.print(o); }

//----------------------- Random integers ----------

#ifdef __GMP_PLUSPLUS__
inline gmp_randclass& Integer::randstate(long unsigned int seed) {
	static gmp_randclass randstate(GMP_RAND_ALG_DEFAULT,seed);
	return static_cast<gmp_randclass&>(randstate);
}

inline void Integer::seeding(long unsigned int s) {
	Integer::randstate().seed(s) ;
}
#endif

inline Integer& Integer::random (Integer& r, long size)
{
#ifdef __GMP_PLUSPLUS__
    mpz_set( (mpz_ptr) &(r.gmp_rep) , ((mpz_class)Integer::randstate().get_z_bits(size)).get_mpz_t() );
#else
    mpz_random((mpz_ptr) &(r.gmp_rep), size);
#endif
    return r;
}


inline Integer Integer::random(int sz)
{
  Integer res;
  return Integer::random(res, sz);
}

inline Integer Integer::nonzerorandom(int sz) {
    Integer r;
    while(isZero(Integer::random(r, sz) )) {};
    return r;
}

inline Integer& Integer::random (Integer& r, const Integer& similar)
{
// J.G.D. 28/04/2004 : z_range thanks to Jack Dubrois
#ifdef __GMP_PLUSPLUS__
    mpz_set( (mpz_ptr) &(r.gmp_rep) , ((mpz_class)Integer::randstate().get_z_range( (mpz_class)( (mpz_ptr) &(similar.gmp_rep) )             )).get_mpz_t() );
#else
    mpz_random((mpz_ptr) &(r.gmp_rep), mpz_size( (mpz_ptr)&(similar.gmp_rep) ) );
#endif
     return r;
}

inline Integer& Integer::nonzerorandom (Integer& r, const Integer& size) {
    while (isZero(Integer::random(r,size))) {};
    return r;
}


inline Integer& Integer::nonzerorandom (Integer& r, long size)
{    while (isZero(Integer::random(r,size))) {};
    return r;
}
