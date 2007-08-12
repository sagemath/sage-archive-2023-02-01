#ifndef _GIVARO_GFQ1_H_
#define _GIVARO_GFQ1_H_

// ==========================================================================
// file: givgfq.h
// Time-stamp: <06 Jun 07 19:00:34 Jean-Guillaume.Dumas@imag.fr>
// (c) Givaro Team
// date: 1999
// version:
// author: Jean-Guillaume.Dumas
// Description:
//   Arithmetic on GF(p^k), with p a prime number less than 2^15
// ==========================================================================

#include "givaro/givconfig.h"
#include "givaro/givinteger.h"
#include <iostream>
#include <vector>
#include "givaro/giv_randiter.h"
#include "givaro/givpoly1.h"

  // ------------------------------------------------- class GFqDom

template<class TT> class GFqDom {
protected:
    typedef typename Signed_Trait<TT>::unsigned_type UTT;
    typedef TT Rep;
public:
    typedef GFqDom<TT> Self_t;


    typedef Rep Element;
//     class Element {
//     public:
//         mutable Rep _value;

//         Element() {}
//     };

    typedef UTT Residu_t;

        // ----- Representation of vector of the Element
    typedef Rep* Array;
    typedef const Rep* constArray;

        // ----- Constantes
//    const Rep zero;
//    const Rep one;
    Rep zero;
    Rep one;

    GFqDom(): zero(0), _log2pol(0), _pol2log(0),_plus1(0) {}

    GFqDom( const UTT P, const UTT e = 1);

    GFqDom( const UTT P, const UTT e, const std::vector<UTT>& modPoly);

    GFqDom( const GFqDom<TT>& F)
      {
	zero = F.zero;
	one = F.one;
	_characteristic = F._characteristic;
	_dcharacteristic = F._dcharacteristic;
	_inversecharacteristic = F._inversecharacteristic;
	_exponent = F._exponent;
	_q = F._q;
	_qm1 = F._qm1;
	_qm1o2 = F._qm1o2;
	_log2pol = F._log2pol;
	_pol2log = F._pol2log;
	_plus1 = F._plus1;
      }

        // Allows to choose the randomization
        // and therefore the field generator
//     template<class RandIter >
//     GFqDom(RandIter& g, const UTT P, const UTT e = 1);
    ~GFqDom() {};


    GFqDom<TT> operator=( const GFqDom<TT>& F)
      {
	this->zero = F.zero;
	this->one = F.one;
	this->_characteristic = F._characteristic;
	this->_dcharacteristic = F._dcharacteristic;
	this->_inversecharacteristic = F._inversecharacteristic;
	this->_exponent = F._exponent;
	this->_q = F._q;
	this->_qm1 = F._qm1;
	this->_qm1o2 = F._qm1o2;
	this->_log2pol = F._log2pol;
	this->_pol2log = F._pol2log;
	this->_plus1 = F._plus1;
	return *this;
      }



        // Access to the modulus, characteristic, size, exponent

    UTT residu() const;
    UTT characteristic() const;
    Integer& characteristic(Integer& p) const{return p=characteristic();}
    UTT generator() const;
    UTT sage_generator() const;
    UTT cardinality() const;
    UTT exponent() const;
    UTT size() const;

        // Initialization of Elements
    Rep& init( Rep&) const;
    Rep& init( Rep&, const int) const ;
    Rep& init( Rep&, const unsigned int) const ;
    Rep& init( Rep&, const long) const ;
    Rep& init( Rep&, const unsigned long) const ;
    Rep& init( Rep&, const Integer) const;
    Rep& init( Rep&, const float) const ;
    Rep& init( Rep&, const double) const ;
#ifndef __GIVARO__DONOTUSE_longlong__
    Rep& init( Rep&, const long long) const;
    Rep& init( Rep&, const unsigned long long) const ;
#endif
    Rep& init( Rep& a, std::istream& s ) const { return read(a,s); }
        // Initialization of a polynomial
    template<typename val_t, template<typename V> class Polynomial>
    Rep& init( Rep&, const Polynomial<val_t>&);


        // -- Misc: r <- a mod p
    Rep& assign (Rep&, const Integer) const;
    Rep& assign (Rep&, const Rep) const;
    void assign ( const size_t sz, Array r, constArray a ) const;

        // --- IO methods for the Domain
    std::istream& read ( std::istream& s );
    std::ostream& write( std::ostream& s ) const;
        // --- IO methods for the Elements
    std::istream& read ( std::istream& s, Rep& a ) const;
    std::ostream& write( std::ostream& s, const Rep a ) const;

        // Conversions of the elements
    std::ostream& 	convert(std::ostream& s, const Rep a ) const { return write(s,a); }
    TT 			convert(const Rep) const ;
    long& 		convert(long&, const Rep) const ;
    unsigned long& 	convert(unsigned long&, const Rep) const ;
    int& 		convert(int&, const Rep) const ;
    float& 		convert(float&, const Rep) const ;
    double& 		convert(double&, const Rep) const ;
    unsigned int& 	convert(unsigned int&, const Rep) const ;
    Integer& 		convert(Integer&, const Rep) const ;
#ifndef __GIVARO__DONOTUSE_longlong__
    long long& 		convert(long long&, const Rep) const ;
    unsigned long long& convert(unsigned long long&, const Rep) const ;
#endif

        // Test operators
    inline int operator== (const GFqDom<TT>& a) const;
    inline int operator!= (const GFqDom<TT>& a) const;

        // Miscellaneous functions
    bool areEqual( const Rep&, const Rep&  ) const;
    bool areNEqual ( const Rep , const Rep ) const;
    bool isZero( const Rep ) const;
    bool isnzero( const Rep ) const;
    bool isOne ( const Rep ) const;
    bool isunit ( const Rep ) const; // Element belongs to prime subfield
    size_t length ( const Rep ) const;



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

    Rep& axpy (Rep& r, const Rep a, const Rep b, const Rep c) const;
    void axpy (const size_t sz, Array r, Rep a, constArray x, constArray y) const;
    void axpy (const size_t sz, Array r, Rep a, constArray x, Rep c) const;

        // -- axpyin: r <- r + a * x mod p
    Rep& axpyin (Rep& r, const Rep a, const Rep b) const;
    void axpyin (const size_t sz, Array r, Rep a, constArray x) const;

        // -- axmy: r <- a * x - y mod p
    Rep& axmy (Rep& r, const Rep a, const Rep b, const Rep c) const;
    void axmy (const size_t sz, Array r, Rep a, constArray x, constArray y) const;
    void axmy (const size_t sz, Array r, Rep a, constArray x, Rep c) const;

        // -- amxy: r <- c - a * b mod p
    Rep& amxy (Rep& r, const Rep a, const Rep b, const Rep c) const;

        // -- axmyin: r <- r - a * b mod p
    Rep& axmyin (Rep& r, const Rep a, const Rep b) const;
    void axmyin (const size_t sz, Array r, Rep a, constArray x) const;

//   // -- sqpyin: r <- r + a * a mod p
//     Rep& sqpyin (Rep& r, const Rep a) const;


        // -- axpyin: r <- r - a * b mod p
    Rep& amxyin (Rep& r, const Rep a, const Rep b) const;

        // <- \sum_i a[i], return 1 if a.size() ==0,
    void reduceadd ( Rep& r, const size_t sz, constArray a ) const;

        // <- \prod_i a[i], return 1 if a.size() ==0,
    void reducemul ( Rep& r, const size_t sz, constArray a ) const;

        // <- \sum_i a[i] * b[i]
    Rep& dotprod ( Rep& r, const size_t sz, constArray a, constArray b ) const;

        // ----- random generators
   // ----- random generators

    template<class RandIter> Rep& random(RandIter& g, Rep& r) const ;
    template<class RandIter> Rep& random(RandIter& g, Rep& r, long s) const ;
    template<class RandIter> Rep& random(RandIter& g, Rep& r, const Rep& b) const ;
    template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const ;
    template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, long s) const ;
    template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const ;

    typedef GIV_randIter< GFqDom<TT> , Rep> randIter;

//         // - Set to a non zero random value
//     void set_nrand(Rep&) const;
//         // - Set to a random value
//     void set_rand(Rep&) const;

#ifdef __GIVARO_COUNT__
    void clear() {
        _add_count = 0;
        _mul_count = 0;
        _neg_count = 0;
        _div_count = 0;
        _sub_count = 0;
        _inv_count = 0;
        _add_call = 0;
        _mul_call = 0;
        _neg_call = 0;
        _div_call = 0;
        _sub_call = 0;
        _inv_call = 0;
    }

    void info() const {
        cerr << "Mul Call: " << _mul_call << ", real: " << _mul_count << endl;
        cerr << "Add Call: " << _add_call << ", real: " << _add_count << endl;
        cerr << "Div Call: " << _div_call << ", real: " << _div_count << endl;
        cerr << "Sub Call: " << _sub_call << ", real: " << _sub_count << endl;
        cerr << "Neg Call: " << _neg_call << ", real: " << _neg_count << endl;
        cerr << "Inv Call: " << _inv_call << ", real: " << _inv_count << endl;
    }
#endif


protected:
//    const UTT _characteristic;
//    const UTT _exponent;
//    const UTT _q;
//    const UTT _qm1;
//    const UTT _qm1o2;
    UTT _characteristic;
    UTT _exponent;
    UTT _q;
    UTT _qm1;
    UTT _qm1o2;

        // G is a generator of GF(q)
        // p is GF(q)'s characteristic
        // log2pol[ i ] = G^i(p)
        // pol2log[ j ] = i such that log2pol[i] = j
        // plus1[i] = k such that G^i + 1 = G^k
    std::vector<UTT> _log2pol;
    std::vector<UTT> _pol2log;
    std::vector<TT> _plus1;

 private:
    double _dcharacteristic;
    double _inversecharacteristic;

#ifdef __GIVARO_COUNT__
static    long long _add_count;
static    long long _mul_count;
static    long long _neg_count;
static    long long _div_count;
static    long long _sub_count;
static    long long _inv_count;
static    long long _add_call;
static    long long _mul_call;
static    long long _neg_call;
static    long long _div_call;
static    long long _sub_call;
static    long long _inv_call;
#endif

    static void Init();
    static void End();
};


#include "givaro/givgfq.inl"

#endif
