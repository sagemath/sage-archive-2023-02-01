#ifndef _GIVARO_ZPZ16STD_H_
#define _GIVARO_ZPZ16STD_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16std.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givzpz16std.h,v 1.9 2006-07-21 08:03:26 jgdumas Exp $
// ==========================================================================
//
//  Modified by Pascal Giorgi on 2002/02/13  (pascal.giorgi@ens-lyon.fr)
//
// Description:
//   Arithmetic on Z/pZ, with p a prime number less than 2^16
//   Modulo typedef is a signed long number. In case it was modified
//   then bezout algorithm must be changed (coefficient can be negative).


#ifdef MACOSX
#  include <sys/types.h> // needed on MacOS X 10.5 for uint type
#endif

/* Thierry -> JG: Constantes necessaires:
   (typedef) uint16 : type des int sur 16bits non signe
   (typedef) uint32 : type des int sur 16bits non signe
   (typedef) int16 : type des int sur 16bits signe
   (typedef) int32 : type des int sur 16bits signe
   (#define) GIVARO_MAXUINT16: 2^16 -1
   (#define) GIVARO_MAXUINT32: 2^32 -1
   (#define) GIVARO_MAXULONG: val max en unsigned long
*/
/*
#define GIVARO_BITS_PER_LONGINT        32
#define GIVARO_BITS_PER_INT            32
#define GIVARO_BITS_PER_SHORTINT       16
#define GIVARO_BITS_PER_CHAR           16
typedef char    int8;
typedef short   int16;
typedef int     int32;
typedef unsigned char   uint8;
typedef unsigned short  uint16;
typedef unsigned int    uint32;

#define GIVARO_MAXUINT8                255U            // 2^8-1
#define GIVARO_MAXUINT16               65535U          // 2^16-1
#define GIVARO_MAXUINT32               4294967295U     // 2^32-1
#define GIVARO_MAXULONG                4294967295U     // 2^32-1
*/
#include "givaro/givbasictype.h"

/*
   Classes d'erreurs:
   GivError::throw_error + des classes d'exception dont les cstos prennent des chaines:
     * GivMathDivZero( " ... " )
     * GivBadFormat( " ... " )
*/
#include "givaro/giverror.h"
#include "givaro/givzpz.h"
#include "givaro/giv_randiter.h"

// ==========================================================================
// -- This class implement the standard arithmetic with Modulo Elements:
// - The representation of an integer a in Zpz is the value a % p
// ==========================================================================

template<>
class ZpzDom<Std16> {
public:
  // ----- Exported Types and constantes
  typedef uint16 Residu_t;                    // - type to store residue
  enum { size_rep = sizeof(Residu_t) };      // - size of the storage type
  // ----- Representation of Element of the domain ZpzDom
  typedef int16 Rep;
  typedef int16 Element;

  // ----- Representation of vector of the Element
  typedef Rep* Array;
  typedef const Rep* constArray;

  // ----- Constantes
  const Rep zero;
  const Rep one;

  // ----- Constructor
  ZpzDom() : zero(0), one(1), _p(0) {}
  ZpzDom( Residu_t p ) : zero(0), one(1), _p(p) {}
  ZpzDom( const ZpzDom<Std16>& F) : zero(0), one(1), _p(F._p) {}


  int operator==( const ZpzDom<Std16>& BC) const { return _p == BC._p;}
  int operator!=( const ZpzDom<Std16>& BC) const { return _p != BC._p;}

  ZpzDom<Std16>& operator=( const ZpzDom<Std16>& F) { this->_p = F._p; return *this;}

  // ----- Access to the modulus
  Residu_t residu() const;
  Residu_t size() const { return _p; }
  Residu_t characteristic() const { return _p; }
  Integer& characteristic( Integer& p) const { return p=_p; }
  Rep access( const Rep a ) const { return a; }

  // ----- Convert from Element to int
    int16& convert( int16& x , const Rep a) const { return x=(int16)(a);}
    uint16& convert( uint16& x , const Rep a) const { return x=(uint16)(a);}
    unsigned long& convert( unsigned long& x , const Rep a) const { return x=(unsigned long)(a);}
    double& convert( double& x , const Rep a) const { return x=(double)(a);}
    int& convert( int& x , const Rep a) const { return x=int(a);}
    Integer& convert(Integer& i, const Rep a) const {
        unsigned long ur;
        return i = (Integer)convert(ur, a);
    }



  // ----- Access to the modulus
  Rep& init( Rep& a ) const;
  void init( const size_t, Array a, constArray b ) const;
  Rep& init( Rep& a, const long i) const ;
  Rep& init( Rep& a, const unsigned long i) const ;
  Rep& init( Rep& a, const int i) const ;
  Rep& init( Rep& a, const unsigned int i) const ;
  Rep& init( Rep& a, const double i) const ;
  Rep& init( Rep& a, const float i) const ;
  Rep& init( Rep& a, const Integer& i) const ;

  // ----- Misc methods
  int areEqual( const  Rep, const Rep) const;
  int areNEqual( const Rep, const Rep) const;
  int isZero( const Rep a ) const;
  int isnzero( const Rep a ) const;
  int isOne ( const Rep a ) const;
  size_t length ( const Rep a ) const;

  // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
  Rep& mul (Rep& r, const Rep a, const Rep b) const;
  Rep& div (Rep& r, const Rep a, const Rep b) const;
  Rep& add (Rep& r, const Rep a, const Rep b) const;
  Rep& sub (Rep& r, const Rep a, const Rep b) const;
  Rep& neg (Rep& r, const Rep a) const;
  Rep& inv (Rep& r, const Rep a) const;

  Rep& mulin (Rep& r, const Rep a) const;
  Rep& divin (Rep& r, const Rep a) const;
  Rep& addin (Rep& r, const Rep a) const;
  Rep& subin (Rep& r, const Rep a) const;
  Rep& negin (Rep& r) const;
  Rep& invin (Rep& r) const;

  // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
  void mul (const size_t sz, Array r, constArray a, constArray b) const;
  void mul (const size_t sz, Array r, constArray a, Rep b) const;

  void div (const size_t sz, Array r, constArray a, constArray b) const;
  void div (const size_t sz, Array r, constArray a, Rep b) const;

  void add (const size_t sz, Array r, constArray a, constArray b) const;
  void add (const size_t sz, Array r, constArray a, Rep b) const;

  void sub (const size_t sz, Array r, constArray a, constArray b) const;
  void sub (const size_t sz, Array r, constArray a, Rep b) const;

  void neg (const size_t sz, Array r, constArray a) const;
  void inv (const size_t sz, Array r, constArray a) const;

  // -- axpy: r <- a * x + y mod p
  Rep& axpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
  void axpy
   (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
  // -- axpyin: r <- r + a * x mod p
  Rep& axpyin(Rep& r, const Rep a, const Rep b) const;
  void axpyin
   (const size_t sz, Array r, constArray a, constArray x) const;

  // -- amxy: r <- c - a * b mod p
  Rep& amxy (Rep& r, const Rep a, const Rep b, const Rep c) const;

  // -- axmy: r <- a * x - y mod p
  Rep& axmy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
  void axmy
   (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
  // -- axmyin: r <- r - a * x mod p
  Rep& axmyin(Rep& r, const Rep a, const Rep b) const;
  void axmyin
   (const size_t sz, Array r, constArray a, constArray x) const;

  // -- Misc: r <- a mod p
  void assign ( const size_t sz, Array r, constArray a ) const;
/* JGD 26.10.99
  void assign ( Rep& r, const Rep a) const;
  void assign ( Rep& r, const long a ) const;
  void assign ( Rep& r, const unsigned long a ) const;
  void assign ( Rep& r, const int a ) const;
  void assign ( Rep& r, const unsigned int a ) const;
*/
  Rep& assign ( Rep& r, const Rep a) const;
  Rep& assign ( Rep& r, const long a ) const;
  Rep& assign ( Rep& r, const unsigned long a ) const;
  Rep& assign ( Rep& r, const int a ) const;
  Rep& assign ( Rep& r, const unsigned int a ) const;
   // ----- random generators
    template< class RandIter > Rep& random(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, const Rep& b) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, const Rep& b) const ;

    typedef GIV_randIter< ZpzDom<Std16>, Rep > randIter;

  // <- \sum_i a[i], return 1 if a.size() ==0,
  void reduceadd ( Rep& r, const size_t sz, constArray a ) const;

  // <- \prod_i a[i], return 1 if a.size() ==0,
  void reducemul ( Rep& r, const size_t sz, constArray a ) const;

  // <- \sum_i a[i] * b[i]
  void dotprod ( Rep& r, const size_t sz, constArray a, constArray b ) const;
  void dotprod ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const;

  // ----- a -> r: uint16 to double
  void i2d ( const size_t sz, double* r, constArray a ) const;

  // ----- a -> r % p: double to uint16 % p
  void d2i ( const size_t sz, Array r, const double* a ) const;

  // --- IO methods
  std::istream& read ( std::istream& s );
  std::ostream& write( std::ostream& s ) const;
  std::istream& read ( std::istream& s, Rep& a ) const;
  std::ostream& write( std::ostream& s, const Rep a ) const;

protected:
  // -- based for modular inverse, d = a*u + b*v
//   static const int32 gcdext ( int32& u, int32& v, const int32 a, const int32 b );
  int32& gcdext (int32& d, int32& u, int32& v, const int32 a, const int32 b ) const;
  int32& invext (int32& u, const int32 a, const int32 b ) const;

protected:
  // -- data representation of the domain:
  Residu_t _p;

  static void Init();
  static void End();
};


#include "givaro/givzpz16std.inl"

#endif
