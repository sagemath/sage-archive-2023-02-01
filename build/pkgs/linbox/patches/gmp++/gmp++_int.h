#ifndef _GMPplusplus_INTEGER_H_
#define _GMPplusplus_INTEGER_H_
// ========================================================================
// Copyright(c)'2001-2007 by LinBox Team
// see the copyright file.
// Authors: M. Samama, T. Gautier, JG. Dumas
// Time-stamp: <11 Jan 07 14:52:53 Jean-Guillaume.Dumas@imag.fr>
// ========================================================================
// Description:
// Integer class definition based on Gmp (>V2.0 or 1.3.2)

// Core gmp++_int.h
#include <vector>
#include <list>
#include <string>


#ifdef __USE_64_bits__
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif

#ifdef __USE_ISOC99
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif


//------------------------------------------------------ Friend Integer
class Integer;

int 		compare(const Integer& a, const Integer& b);
int 		absCompare(const Integer& a, const Integer& b);
Integer& 	inv (Integer& u, const Integer& a, const Integer& b);
Integer 	gcd (const Integer& a, const Integer& b);
Integer 	gcd (const Integer& a, const Integer& b, Integer& u, Integer& v);
Integer& 	gcd (Integer& g, const Integer& a, const Integer& b);
Integer& 	gcd (Integer& g, const Integer& a, const Integer& b, Integer& u, Integer& v);
Integer 	pp( const Integer& P, const Integer& Q );
Integer& 	lcm (Integer& g, const Integer& a, const Integer& b);
Integer 	lcm (const Integer& a, const Integer& b);
Integer& 	pow(Integer& Res, const Integer& n, const long l);
Integer& 	pow(Integer& Res, const unsigned long n, const unsigned long l);
Integer& 	pow(Integer& Res, const Integer& n, const unsigned long l);
Integer& 	pow(Integer& Res, const Integer& n, const int l) ;
Integer& 	pow(Integer& Res, const Integer& n, const unsigned int l) ;
Integer 	pow(const Integer& n, const long l);
Integer 	pow(const Integer& n, const unsigned long l);
Integer 	pow(const Integer& n, const int l) ;
Integer 	pow(const Integer& n, const unsigned int l);
Integer& 	powmod(Integer& Res, const Integer& n, const unsigned long e, const Integer& m);
Integer& 	powmod(Integer& Res, const Integer& n, const long e, const Integer& m);
Integer& 	powmod(Integer& Res, const Integer& n, const unsigned int e, const Integer& m) ;
Integer& 	powmod(Integer& Res, const Integer& n, const int e, const Integer& m)  ;
Integer& 	powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);
Integer 	powmod(const Integer& n, const unsigned long e, const Integer& m);
Integer 	powmod(const Integer& n, const long e, const Integer& m);
Integer 	powmod(const Integer& n, const unsigned int e, const Integer& m) ;
Integer 	powmod(const Integer& n, const int e, const Integer& m) ;
Integer 	powmod(const Integer& n, const Integer& e, const Integer& m);
Integer 	fact ( unsigned long l);
Integer 	sqrt(const Integer& p);
Integer 	sqrtrem(const Integer& p, Integer& rem);
Integer& 	sqrt(Integer& r, const Integer& p);
Integer& 	sqrtrem(Integer& r, const Integer& p, Integer& rem);
bool 		root(Integer& q, const Integer&, unsigned int n);
long 		logp(const Integer& a, const Integer& p) ;
double 		logtwo(const Integer& a) ;
void 		swap(Integer& , Integer&);
int 		sign   (const Integer& a);
int 		isZero (const Integer& a);
int 		isOne  (const Integer& a);
int 		isperfectpower  (const Integer& );
Integer 	abs(const Integer& n);
Integer& 	prevprime(Integer&, const Integer& p);
Integer& 	nextprime(Integer&, const Integer& p);
int 		probab_prime(const Integer& p);
int 		probab_prime(const Integer& p, int r);
int 		jacobi(const Integer& u, const Integer& v) ;
int 		legendre(const Integer& u, const Integer& v) ;
unsigned long 	length (const Integer& a);
std::istream& 	operator >> (std::istream &i, Integer& n);
std::ostream& 	operator << (std::ostream &o, const Integer& n);
std::ostream& 	absOutput (std::ostream &o, const Integer& n);
void 		importWords(Integer&, size_t, int, int, int, size_t, const void*);


//------------------------------------------------------ Class Integer
class Integer {

public:
  typedef std::vector<mp_limb_t> vect_t;
  Integer( const std::vector<mp_limb_t>& vect_t );
  //--------------------------------------cstors & dstors
  Integer(int n = 0);
  Integer(long n);
  Integer(unsigned char n);
  Integer(unsigned int n);
  Integer(unsigned long n);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
  Integer(long long n);
  Integer(unsigned long long n);
#endif
  Integer(double d);
  Integer(const char *s);
  Integer(const Integer& n);
  ~Integer();

  //------------------------------------- predefined null and one
  static const Integer zero;
  static const Integer one;

  // -- Assignment and copy operators
  Integer& operator = (const Integer& n);
  Integer& logcpy(const Integer& n);
  Integer& copy(const Integer& n);

  //------------------Equalities and inequalities between integers and longs
  int operator != (const int l) const;
  int operator != (const long l) const;
  int operator != (const unsigned long l) const;

  friend int compare(const Integer& a, const Integer& b);
  friend int absCompare(const Integer& a, const Integer& b);

  int operator > (const int l) const;
  int operator > (const long l) const;
  int operator > (const unsigned long l) const;
  int operator < (const int l) const;
  int operator < (const long l) const;
  int operator < (const unsigned long l) const;

  //------------------ Bit logic
  Integer operator^ (const Integer&);   // XOR
  Integer operator| (const Integer&);   // OR
  Integer operator& (const Integer&);   // AND
  Integer operator ~ () const;   // 1 complement
  Integer& operator^= (const Integer&);   // XOR
  Integer& operator|= (const Integer&);   // OR
  Integer& operator&= (const Integer&);   // AND
  Integer operator<< (int l) const; // lshift
  Integer operator>> (int l) const; // rshift
  Integer operator<< (long l) const; // lshift
  Integer operator>> (long l) const; // rshift
  Integer operator<< (unsigned int l) const; // lshift
  Integer operator>> (unsigned int l) const; // rshift
  Integer operator<< (unsigned long l) const; // lshift
  Integer operator>> (unsigned long l) const; // rshift
  Integer& operator<<= (int l) ; // lshift
  Integer& operator>>= (int l) ; // rshift
  Integer& operator<<= (long l) ; // lshift
  Integer& operator>>= (long l) ; // rshift
  Integer& operator<<= (unsigned int l) ; // lshift
  Integer& operator>>= (unsigned int l) ; // rshift
  Integer& operator<<= (unsigned long l) ; // lshift
  Integer& operator>>= (unsigned long l) ; // rshift

  //----------------Elementary arithmetic between Integers & longs
  Integer  operator + (const Integer& n) const;
  Integer  operator + (const unsigned long l) const;
  Integer  operator + (const long l) const;
  Integer& operator += (const Integer& n);
  Integer& operator += (const unsigned long l);
  Integer& operator += (const long l);
  template<class XXX> Integer& operator +=(const XXX& x) { return this->operator += ( (Integer)x ); }

  Integer  operator - (const Integer& n) const;
  Integer  operator - (const unsigned long l) const;
  Integer  operator - (const long l) const;
  Integer& operator -= (const Integer& n);
  Integer& operator -= (const unsigned long l);
  Integer& operator -= (const long l);
  template<class XXX> Integer& operator -=(const XXX& x) { return this->operator -= ( (Integer)x ); }
  Integer  operator -() const;

  Integer  operator * (const Integer& n) const;
  Integer  operator * (const unsigned long l) const;
  Integer  operator * (const long l) const;
  Integer& operator *= (const Integer& n);
  Integer& operator *= (const unsigned long l);
  Integer& operator *= (const long l);
  template<class XXX> Integer& operator *=(const XXX& x) { return this->operator *= ( (Integer)x ); }

  // -- Euclidian division of a/b: returns q or r such that
  // - a=b*q + r, with |r| < |b|, a*r >=0
  Integer  operator /  (const Integer& n) const;
  Integer  operator /  (const unsigned long l) const;
  Integer  operator /  (const long l) const;
  Integer& operator /= (const Integer& n);
  Integer& operator /= (const unsigned long l);
  Integer& operator /= (const long l);
  template<class XXX> Integer& operator /=(const XXX& x) { return this->operator /= ( (Integer)x ); }

  Integer  operator % (const Integer& n) const;
  unsigned long  operator % (const unsigned long l) const;
  long  operator % (const long l) const;
  unsigned short  operator % (const unsigned short l) const { return (unsigned short) ( this->operator % ( (unsigned long)l ) ); }
  template<class XXX> XXX operator %(const XXX& x) const { return (XXX)this->operator % ( Integer(x) ); }
  Integer& operator %= (const Integer& n);
  Integer& operator %= (const unsigned long l);
  Integer& operator %= (const long l);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
  Integer& operator %= (const long long l) { return *this %= (Integer)l; }
  Integer& operator %= (const unsigned long long l) { return *this %= (Integer)l; }
  long long operator % (const long long l) const;
  unsigned long long operator % (const unsigned long long l) const;
#endif
  template<class XXX> Integer& operator %=(const XXX& x) { return this->operator %= ( (Integer)x ); }

  // - Methods
static Integer& addin (Integer& res, const Integer& n);
static Integer& addin (Integer& res, const long n);
static Integer& addin (Integer& res, const unsigned long n);
static Integer& add   (Integer& res, const Integer& n1, const Integer& n2);
static Integer& add   (Integer& res, const Integer& n1, const long n2);
static Integer& add   (Integer& res, const Integer& n1, const unsigned long n2);

static Integer& subin (Integer& res, const Integer& n);
static Integer& subin (Integer& res, const long n);
static Integer& subin (Integer& res, const unsigned long n);
static Integer& sub   (Integer& res, const Integer& n1, const Integer& n2);
static Integer& sub   (Integer& res, const Integer& n1, const long n2);
static Integer& sub   (Integer& res, const Integer& n1, const unsigned long n2);
static Integer& negin (Integer& res);
static Integer& neg   (Integer& res, const Integer& n);

static Integer& mulin (Integer& res, const Integer& n);
static Integer& mulin (Integer& res, const long n);
static Integer& mulin (Integer& res, const unsigned long n);
static Integer& mul   (Integer& res, const Integer& n1, const Integer& n2);
static Integer& mul   (Integer& res, const Integer& n1, const long n2);
static Integer& mul   (Integer& res, const Integer& n1, const unsigned long n2);
static Integer& axpy   (Integer& res, const Integer& a, const Integer& x, const Integer& y );
static Integer& axpyin   (Integer& res, const Integer& a, const Integer& x);
static Integer& amxy   (Integer& res, const Integer& a, const Integer& x, const Integer& y );
static Integer& axmy   (Integer& res, const Integer& a, const Integer& x, const Integer& y );
static Integer& axmyin   (Integer& res, const Integer& a, const Integer& x);

static Integer& divin (Integer& q, const Integer& n);
static Integer& divin (Integer& q, const long n);
static Integer& divin (Integer& q, const unsigned long n);
static Integer& div   (Integer& q, const Integer& n1, const Integer& n2);
static Integer& div   (Integer& q, const Integer& n1, const long n2);
static Integer& div   (Integer& q, const Integer& n1, const unsigned long n2);
static Integer& divexact  (Integer& q, const Integer& n1, const Integer& n2);
static Integer  divexact  (const Integer& n1, const Integer& n2);

static Integer& modin (Integer& r, const Integer& n);
static Integer& modin (Integer& r, const long n);
static Integer& modin (Integer& r, const unsigned long n);
static Integer& mod   (Integer& r, const Integer& n1, const Integer& n2);
static Integer& mod   (Integer& r, const Integer& n1, const long n2);
static Integer& mod   (Integer& r, const Integer& n1, const unsigned long n2);

  // -- return q, the quotient
static Integer& divmod   (Integer& q, Integer& r, const Integer& n1, const Integer& n2);
static Integer& divmod   (Integer& q, long& r, const Integer& n1, const long n2);
static Integer& divmod   (Integer& q, unsigned long& r, const Integer& n1, const unsigned long n2);


  //------------------------------------- Arithmetic functions
  friend Integer& inv (Integer& u, const Integer& a, const Integer& b);
  friend Integer gcd (const Integer& a, const Integer& b);
  friend Integer gcd (const Integer& a, const Integer& b,
                            Integer& u, Integer& v);
  friend Integer& gcd (Integer& g, const Integer& a, const Integer& b);
  friend Integer& gcd (Integer& g, const Integer& a, const Integer& b,
                            Integer& u, Integer& v);

  friend Integer pp( const Integer& P, const Integer& Q );

  friend Integer& lcm (Integer& g, const Integer& a, const Integer& b);
  friend Integer lcm (const Integer& a, const Integer& b);

  // - return n^l
  friend Integer& pow(Integer& Res, const Integer& n, const long l);
  friend Integer& pow(Integer& Res, const unsigned long n, const unsigned long l);
  friend Integer& pow(Integer& Res, const Integer& n, const unsigned long l);
  friend Integer& pow(Integer& Res, const Integer& n, const int l) { return pow(Res, n, (long)l ); }
  friend Integer& pow(Integer& Res, const Integer& n, const unsigned int l) { return pow(Res, n, (unsigned long)l ); }
  friend Integer pow(const Integer& n, const long l);
  friend Integer pow(const Integer& n, const unsigned long l);
  friend Integer pow(const Integer& n, const int l) { return pow(n, (long)l ); }
  friend Integer pow(const Integer& n, const unsigned int l) { return pow(n, (unsigned long)l ); }

  // - return n^e % m
  friend Integer& powmod(Integer& Res, const Integer& n, const unsigned long e, const Integer& m);
  friend Integer& powmod(Integer& Res, const Integer& n, const long e, const Integer& m);
  friend Integer& powmod(Integer& Res, const Integer& n, const unsigned int e, const Integer& m) { return powmod(Res, n, (unsigned long)e, m); }
  friend Integer& powmod(Integer& Res, const Integer& n, const int e, const Integer& m)  { return powmod(Res, n, (long)e, m); }
  friend Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);
  friend Integer powmod(const Integer& n, const unsigned long e, const Integer& m);
  friend Integer powmod(const Integer& n, const long e, const Integer& m);
  friend Integer powmod(const Integer& n, const unsigned int e, const Integer& m) { return powmod(n, (unsigned long)e, m); }
  friend Integer powmod(const Integer& n, const int e, const Integer& m)  { return powmod(n, (long)e, m); }
  friend Integer powmod(const Integer& n, const Integer& e, const Integer& m);

  friend Integer fact ( unsigned long l);

  friend Integer sqrt(const Integer& p);
  friend Integer sqrtrem(const Integer& p, Integer& rem);
  friend Integer& sqrt(Integer& r, const Integer& p);
  friend Integer& sqrtrem(Integer& r, const Integer& p, Integer& rem);


  friend bool root(Integer& q, const Integer&, unsigned int n);
  friend long logp(const Integer& a, const Integer& p) ;
  friend double logtwo(const Integer& a) ;

  //-----------------------------------------Miscellaneous
  friend void swap(Integer& , Integer&);

  friend inline int sign   (const Integer& a);
  friend inline int isZero (const Integer& a);
  friend inline int isOne  (const Integer& a);
  friend int isperfectpower  (const Integer& );

  friend Integer abs(const Integer& n);

  friend Integer& prevprime(Integer&, const Integer& p);
  friend Integer& nextprime(Integer&, const Integer& p);
  friend int probab_prime(const Integer& p);
  friend int probab_prime(const Integer& p, int r);
  friend int jacobi(const Integer& u, const Integer& v) ;
  friend int legendre(const Integer& u, const Integer& v) ;

  Integer& operator++() { return *this+=1UL; } // prefix
  Integer& operator--() { return *this-=1UL; } // prefix

  // - return the size in byte
  friend inline unsigned long length (const Integer& a);
  // - return the size in word.
  size_t size() const;
  // - return the size in base B (always exact if B is a power of two)
  size_t size_in_base(int B) const;
  // - return the size in bit.
  size_t bitsize() const;
  // - return the i-th word of the integer. Word 0 is lowest word.
  unsigned long operator[](size_t i) const;

  // -- Convert an Integer to a basic C++ type
  // -- Cast operators
  operator short() const { return (int) *this; }
  operator unsigned short() const { return (unsigned int) *this; }
  operator unsigned char() const { return (unsigned int) *this; }
  operator unsigned int() const ;
  operator int() const ;
  operator signed char() const { return (int) *this; }
  operator unsigned long() const ;
  operator long() const ;
#ifndef __GIVARO__DONOTUSE_longlong__
  operator unsigned long long() const ;
  operator long long() const ;
#endif
  operator std::string() const ;
  operator float() const ;
  operator double() const ;
  operator vect_t() const ;

  //--------------------Random Iterators
  // -- return a random number with sz machine word.
  // -- To be improved.
#ifdef __GMP_PLUSPLUS__
    static void seeding(unsigned long int s=0);
    static gmp_randclass& randstate(unsigned long int s=0);
#endif
    static Integer  random(int sz=1 );
    static Integer  nonzerorandom(int sz=1 );
    static Integer& random(Integer& r, const Integer& size );
    static Integer& nonzerorandom(Integer& r, const Integer& size );
    static Integer& random(Integer& r, long size =1 );
    static Integer& nonzerorandom(Integer& r, long size =1 );
  //----------------------------------------------I/O

  friend std::istream& operator >> (std::istream &i, Integer& n);
  friend std::ostream& operator << (std::ostream &o, const Integer& n);
  friend std::ostream& absOutput (std::ostream &o, const Integer& n);

  friend void importWords(Integer&, size_t, int, int, int, size_t, const void*);

  std::ostream& print( std::ostream& o ) const;

protected:

    typedef MP_INT Rep;

    Rep gmp_rep;

    int priv_sign() const;
    mpz_ptr get_mpz() {return (mpz_ptr)&gmp_rep;}
    const Rep* get_rep() const { return &gmp_rep; }

    // -- Creates a new Integer from a size sz and a array of unsigned long d
    Integer(unsigned long* d, long size);

}; //----------------------------------------------- End of Class Integer



#include "gmp++/gmp++_int.inl"

#endif
