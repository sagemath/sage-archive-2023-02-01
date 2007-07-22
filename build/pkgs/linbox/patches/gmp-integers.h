/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* File: gmp-integers.h
 *  Author: bds
 */

#ifndef __LINBOX_GMP_INTEGERS_H
#define __LINBOX_GMP_INTEGERS_H

#include <linbox/integer.h>
//#include <iostream>
//#include <linbox/util/debug.h>
//#include <linbox/randiter/ntl-ZZ.h>
#include <linbox/field/unparametric.h>
#include <linbox/field/field-traits.h>

namespace LinBox {

	/** wrapper of GMP's integers as a LinBox ring.
	\ingroup ring
	*/
	typedef UnparametricField<integer> GMP_Integers;

	template <class Ring>
    struct ClassifyRing;

	template<>
	struct ClassifyRing<GMP_Integers> {
		typedef RingCategories::IntegerTag categoryTag;
	};

	template <>
	GMP_Integers::Element& GMP_Integers::init(GMP_Integers::Element& x, const integer& y) const {
		return x = y;
	}

#if 0
	class GMP_Integers {

	public:
		typedef NTL_ZZRandIter RandIter;

		typedef NTL::ZZ Element;

		NTL_ZZ(int p = 0, int exp = 1) {
			if( p != 0 ) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be 0 (no modulus)");
			if( exp != 1 ) throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 1");
		}

		inline integer& cardinality (integer& c) const {
			return c = -1;
		}

		inline integer& characteristic (integer& c)const   {
			return c = 0;
		}

		std::ostream& write (std::ostream& out) const  {
			return out << "NTL ZZ Ring";
		}

		std::istream& read (std::istream& in) const  {
			return in;
		}

		/** @brief
		 *  Init x from y.
		 */
		template<class Element2>
		inline Element& init (Element& x,  const Element2& y) const  {

			NTL::conv (x, y);

			return x;
		}

		/** @brief
		 *   Init from a NTL::ZZ
                 */
                inline Element& init (Element& x, const Element& y) const {

			x = y;

			return x;
		}

		/** @brief
		 *   Init from an int64
		 */
		inline Element& init (Element& x, const int64& y) const {
			bool isNeg = false;
			uint64 t;
			if( y < 0 ) {
				isNeg = true;
				t = y * -1;
			}
			else t = y;
			init(x,t);
			if( isNeg ) x *= -1;
			return x;
		}

      /** @brief
      *   Init from a uint64
      */
      inline Element& init (Element& x, const uint64& y) const {
      	uint64 shift = (uint64)1 << 32;
      	uint32 temp = y % shift;
      	NTL::conv (x,temp);
      	x <<= 32;
      	temp = y / shift;
      	x += temp;
      	return x;
      }

		/** @brief
		 *  I don't  know how to init from integer efficiently.
		 */
		 // c_str is safer than data, Z. W and BDS
		inline Element& init (Element& x, const integer& y) const {

			return x=NTL::to_ZZ((std::string(y)).c_str());
		}

		/** @brief
		 *  Convert y to an Element.
		 */
		static inline integer& convert (integer& x, const Element& y){
			bool neg=false;
			if (sign(y) <0)
				neg=true;
			long b = NumBytes(y);
			unsigned char* byteArray;
			byteArray = new unsigned char[(size_t)b ];
			BytesFromZZ(byteArray, y, b);

			integer base(256);
			x= integer(0);

			for(long i = b - 1; i >= 0; --i) {
				x *= base;
				x += integer(byteArray[i]);
			}
			delete [] byteArray;
			if (neg)
				x=-x;
			return x;
		}

		static inline double& convert (double& x, const Element& y){
			return x=NTL::to_double(y);
		}



		/** @brief
		 *  x = y.
		 */
		inline Element&  assign (Element& x, const Element& y)  const {
			return x = y;
		}

		/** @brief
		 *  Test if x == y
		 */
		inline bool areEqual (const Element& x ,const Element& y) const  {
			return x == y;
		}

		/** @brief
		 *  Test if x == 0
		 */
		inline bool isZero (const Element& x) const  {
			return NTL::IsZero (x);
		}

		/** @brief
		 *  Test if x == 1
		 */
		inline bool isOne (const Element& x) const  {
			return NTL::IsOne (x);
		}

		// arithmetic

		/** @brief
		 *  return x = y + z
		 */
		inline Element& add (Element& x, const Element& y, const Element& z) const  {

			NTL::add (x, y, z);

			return x;
		}

		/** @brief
		 *  return x = y - z
		 */
		inline Element& sub (Element& x, const Element& y, const Element& z) const  {

			NTL::sub (x, y, z);

			return x;
		}

		/** @brief
		 *  return x = y * z
		 */
		template <class Int>
		inline Element& mul (Element& x, const Element& y, const Int& z) const  {

			NTL::mul (x, y, z);

			return x;
		}

		/** @brief
		 *  If z divides y, return x = y / z,
		 *  otherwise, throw an exception
		 */
		inline Element& div (Element& x, const Element& y, const Element& z) const {

			Element q, r;

			NTL::DivRem (q, r, y, z);

			if (NTL::IsZero (r))
				return x = q;

			else
				throw PreconditionFailed(__FUNCTION__,__LINE__,"Div: not dividable");
		}

		/** @brief
		 *  If y is a unit, return x = 1 / y,
		 *  otherwise, throw an exception
		 */
		inline Element& inv (Element& x, const Element& y) const {

			if ( NTL::IsOne (y)) return x = y;

			else if ( NTL::IsOne (-y)) return x = y;

			else
				throw PreconditionFailed(__FUNCTION__,__LINE__,"Inv: Not invertible");
		}

		/** @brief
		 *  return x = -y;
		 */
		inline Element& neg (Element& x, const Element& y) const  {

			NTL::negate (x, y);

			return x;
		}


		/** @brief
		 *  return r = a x + y
		 */

		template <class Int>
		inline Element& axpy (Element& r, const Element& a, const Int& x, const Element& y) const  {

			NTL::mul (r, a, x);

			return r += y;
		}


		// inplace operator

		/** @brief
		 *  return x += y;
		 */
		inline Element& addin (Element& x, const Element& y) const {

			return x += y;
		}

		/** @brief
		 *  return x -= y;
		 */
		inline Element& subin (Element& x, const Element& y)  const {

			return x -= y;
		}

		/** @brief
		 *  return x *= y;
		 */
		template<class Int>
		inline Element& mulin (Element& x, const Int& y)  const {

			return x *= y;
		}

		/** @brief
		 *  If y divides x, return x /= y,
		 *  otherwise throw an exception
		 */
		inline Element& divin (Element& x, const Element& y) const {

			div (x, x, y);

			return x;
		}

		/** @brief
		 *  If x is a unit, x = 1 / x,
		 *  otherwise, throw an exception.
		 */
		inline Element& invin (Element& x) {

			if (NTL::IsOne (x)) return x;

			else if (NTL::IsOne (-x)) return x;

			else throw PreconditionFailed(__FUNCTION__,__LINE__,"Div: not dividable");
		}

		/** @brief
		 *  return x = -x;
		 */
		inline Element& negin (Element& x) const  {

			NTL::negate (x, x);

			return x;
		}

		/** @brief
		 *  return r += a x
		 */
		template <class Int>
		inline Element& axpyin (Element& r, const Element& a, const Int& x) const  {

			return r += a * x;
		}


		// IO

		/** @brief
		 *  out << y;
		 */
		std::ostream& write(std::ostream& out,const Element& y) const  {

			out << y;

			return out;
		}


		/** @brief
		 *  read x from istream in
		 */
		std::istream& read(std::istream& in, Element& x) const {

			return in >> x;
		}


		/** some PIR function
		 */

		/** @brief
		 *  Test if x is a unit.
		 */
		inline bool isUnit (const Element& x) const {

			return (NTL::IsOne (x) || NTL::IsOne (-x));
		}

		/** @brief
		 *  return g = gcd (a, b)
		 */
		inline Element& gcd (Element& g, const Element& a, const Element& b) const {

			NTL::GCD (g, a, b);

			return g;
		}

		/** @brief
		 *  return g = gcd (g, b)
		 */
		inline Element& gcdin (Element& g, const Element& b) const {

			NTL::GCD (g, g, b);

			return g;
		}

		/** @brief
		 *  g = gcd(a, b) = a*s + b*t.
		 *  The coefficients s and t are defined according to the standard
		 *  Euclidean algorithm applied to |a| and |b|, with the signs then
		 *  adjusted according to the signs of a and b.
		 */
		inline Element& xgcd (Element& g, Element& s, Element& t, const Element& a, const Element& b)const {

			NTL::XGCD (g,s,t,a,b);

			return g;
		}

		/** @brief
		 *  c = lcm (a, b)
		 */
		inline Element& lcm (Element& c, const Element& a, const Element& b) const {


			if (NTL::IsZero (a) || NTL::IsZero (b)) return c = NTL::ZZ::zero();

			else {
				Element g;

				NTL::GCD (g, a, b);

				NTL::mul (c, a, b);

				c /= g;

				NTL::abs (c, c);

				return c;
			}
		}

		/** @brief
		 *  l = lcm (l, b)
		 */
		inline Element& lcmin (Element& l, const Element& b) const {

			if (NTL::IsZero (l) || NTL::IsZero (b))

				return l = NTL::ZZ::zero();

			else {

				Element g;

				NTL::GCD (g, l, b);

				l *= b;

				l /= g;

				NTL::abs (l, l);

				return l;
			}
		}





		// some specail function

		/** @brief
		 *  x = floor ( sqrt(y)).
		 */

		inline Element& sqrt (Element& x, const Element& y) const  {

			NTL::SqrRoot(x,y);

			return x;
		}

		/** @brief
		 *  Requires 0 <= x < m, m > 2 * a_bound * b_bound,
		 *  a_bound >= 0, b_bound > 0
		 *   This routine either returns 0, leaving a and b unchanged,
		 *   or returns 1 and sets a and b so that
		 *  (1) a = b x (mod m),
		 *  (2) |a| <= a_bound, 0 < b <= b_bound, and
		 *  (3) gcd(m, b) = gcd(a, b).
		 */

		inline long reconstructRational (Element& a, Element& b, const Element& x, const Element& m,
							const Element& a_bound, const Element& b_bound) const {

			return NTL::ReconstructRational(a,b,x,m,a_bound,b_bound);
		}


		/** @brief
		 *  q = floor (x/y);
		 */
		inline Element& quo (Element& q, const Element& a, const Element& b) const {

			NTL::div (q, a, b);

			return q;
		}

		/** @brief
		 *  r = reminder of  a / b
		 */
		inline Element& rem (Element& r, const Element& a, const Element& b) const  {

			NTL::rem (r, a, b);

			return r;
		}

		/** @brief
		 *  a = quotient (a, b)
		 */
		inline Element& quoin (Element& a, const Element& b) const  {

			return a /= b;

		}

		/** @brief
		 *  a = quotient (a, b)
		 */
		inline Element& remin (Element& x, const Element& y)  const {
			return x %= y;
		}


		/** @brief
		 * q = [a/b], r = a - b*q
		 * |r| < |b|, and if r != 0, sign(r) = sign(b)
		 */
		inline void quoRem (Element& q, Element& r, const Element& a, const Element& b) const {

			NTL::DivRem(q,r,a,b);
		}

		/** @brief
		 *  Test if b | a.
		 */
		inline bool isDivisor (const Element& a, const Element& b) const {

			if ( NTL::IsZero (a) ) return true;

			else if (NTL::IsZero (b)) return false;

			else {
				Element r;

				NTL::rem (r, a, b); //weird order changed, dpritcha 2004-07-19

				return NTL::IsZero (r);
			}
		}

		/** compare two elements, a and b
		  * return 1, if a > b
		  * return 0, if a = b;
		  * return -1. if a < b
		  */
		inline long compare (const Element& a, const Element& b) const {

			return NTL::compare (a, b);
		}

		/** return the absolute value
		  * x = abs (a);
		  */
		inline Element& abs (Element& x, const Element& a) const {

			NTL::abs (x, a);

			return x;
		}


		static inline int getMaxModulus() { return 0; } // no modulus

	};


	template<>
	class FieldAXPY<NTL_ZZ>  {
	public:
		typedef NTL_ZZ Field;
		typedef Field::Element Element;

		/** Constructor.
                 * A faxpy object if constructed from a Field and a field element.
                 * Copies of this objects are stored in the faxpy object.
                 * @param F field F in which arithmetic is done
                 */
                FieldAXPY (const Field &F) : _F (F) { _y = 0; }

                /** Copy constructor.
                 * @param faxpy
                 */
                FieldAXPY (const FieldAXPY<Field> &faxpy) : _F (faxpy._F), _y (faxpy._y) {}

                /** Assignment operator
                 * @param faxpy
                 */
                FieldAXPY<Field> &operator = (const FieldAXPY &faxpy)
		{ _y = faxpy._y; return *this; }

                /** Add a*x to y
                 * y += a*x.
                 * @param a constant reference to element a
                 * @param x constant reference to element x
		 * allow optimal multiplication, such as integer * int
                 */
		template<class Element1>
                inline Element& mulacc (const Element &a, const Element1 &x)
		{
			return _y += a * x;
		}

		/** Add a*x to y
                 * y += a*x.
                 * @param a constant reference to element a
                 * @param x constant reference to element x
                 * allow optimal multiplication, such as integer * int
                 */
                template<class Element1>
                inline Element& mulacc (const Element1 &a, const Element &x)
                {
                        return _y += a * x;
                }

		inline Element& mulacc (const Element& a, const Element& b) {

			return _y += a * b;
		}

            inline Element& accumulate (const Element& t) {

			return _y += t;
		}


                /** Retrieve y
                 *
                 * Performs the delayed modding out if necessary
                 */
                inline Element &get (Element &y) { y = _y; return y; }

                /** Assign method.
                 * Stores new field element for arithmetic.
                 * @return reference to self
                 * @param y_init constant reference to element a
                 */
                inline FieldAXPY &assign (const Element& y)
                {
                        _y = y;
                        return *this;
                }

		inline void reset() {
			_y = 0;
		}

	private:

                /// Field in which arithmetic is done
                /// Not sure why it must be mutable, but the compiler complains otherwise
                Field _F;

                /// Field element for arithmetic
                Element _y;

	};

#endif
}

#endif
