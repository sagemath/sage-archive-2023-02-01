/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/modular-float.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2007 Clement Pernet
 * Written by Clement Pernet <cpernet@uwaterloo.ca>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */




#ifndef __LINBOX_MODULAR_FLOAT_H
#define __LINBOX_MODULAR_FLOAT_H


#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <math.h>
#include <linbox/field/field-traits.h>




// Namespace in which all LinBox code resides
namespace LinBox {

	template< class Element >
	class Modular;
	template< class Element >
	class ModularRandIter;

	template< class Field, class RandIter >
	class NonzeroRandIter;

	template <class Ring>
	struct ClassifyRing;
	template <class Element>
	struct ClassifyRing<Modular<Element> >;
	template <>
	struct ClassifyRing<Modular<float> >{
		typedef RingCategories::ModularTag categoryTag;
	};

	class MultiModFloat;

	/// \ingroup field
	template <>
	class Modular<float> : public FieldInterface {

	protected:

		float  modulus;
		unsigned long   lmodulus;

		//float inv_modulus;

	public:
		friend class FieldAXPY<Modular<float> >;
		friend class DotProductDomain<Modular<float> >;
		friend class MultiModFloat;

		typedef float Element;
		typedef ModularRandIter<float> RandIter;
		typedef NonzeroRandIter<Modular<double>, ModularRandIter<double> > NonZeroRandIter;

		static ClassifyRing<Modular<float> >::categoryTag getCategory() {return ClassifyRing<Modular<float> >::categoryTag();}



		Modular () {}

		Modular (int32 p, int exp = 1)  : modulus((float)p), lmodulus(p)//, inv_modulus(1./(float)p)
		{
			if(modulus <= 1)
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			if( exp != 1 ) throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 1");
			integer max;
			if(modulus > (float) FieldTraits<Modular<float> >::maxModulus(max))
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");

		}

		Modular (float p) : modulus(p), lmodulus((unsigned long)p) {
			if( modulus <= 1 )
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			integer max;
			if( modulus > (float) FieldTraits<Modular<float> >::maxModulus(max))
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
		}

		Modular (long int p) :modulus((float)p), lmodulus(p) {
			if( (float) modulus <= 1 )
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			integer max;
			if( (float) modulus > (float) FieldTraits<Modular<float> >::maxModulus(max))
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
		}

		Modular (const integer& p) : modulus((float) p), lmodulus(p) //, inv_modulus(1./(float)p)
		{
			if(modulus <= 1)
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			integer max;
			if(modulus > (float) FieldTraits<Modular<float> >::maxModulus(max))
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");

		}

		Modular(const Modular<float>& mf) : modulus(mf.modulus), lmodulus(mf.lmodulus)//,inv_modulus(mf.inv_modulus)
		{}

		const Modular &operator=(const Modular<float> &F) {
			modulus = F.modulus;
			lmodulus= F.lmodulus;
			//inv_modulus = F.inv_modulus;
			return *this;
		}


		integer &cardinality (integer &c) const{
			return c = integer(modulus);
		}

		integer &characteristic (integer &c) const {
			return c = integer(modulus);
		}

		integer &convert (integer &x, const Element &y) const {
			return x = integer(y);
		}

		float &convert (float &x, const Element& y) const {
			return x=y;
		}

		std::ostream &write (std::ostream &os) const {
			return os << "float mod " << (int)modulus;
		}

		std::istream &read (std::istream &is) {
			is >> modulus;
			if(modulus <= 1)
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
		 	if(modulus > 94906265)
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");

			return is;
		}

		std::ostream &write (std::ostream &os, const Element &x) const {
			return os << x;
		}

		std::istream &read (std::istream &is, Element &x) const {
			integer tmp;
                            // JGD : should'nt it be float tmp ???
			is >> tmp;
			init(x,tmp);
			return is;
		}


		Element &init (Element &x, const integer &y) const  {
// 			return x = (Element)mpz_fdiv_ui(y.get_mpz(),lmodulus );
			return x = (Element)(y%lmodulus);
		}

		inline Element& init(Element& x, float y =0.0) const {

			//float tmp = y;

			/*
			int sign=0;
			if (tmp < 0.0) {
				tmp=-tmp;
				sign=1;
			}
			*/

			//			tmp = floor (y + 0.5);

			//Some odds donot support it. It is in C99.
			//tmp = round (y);

			x = fmod (y, modulus);

			/*
			if (tmp > modulus)
				tmp -= (modulus * floor( tmp*inv_modulus));

			if ( (!tmp) || (tmp == modulus) ){
				return x = 0.0;

			}
			else
				if (sign)
					return x = modulus-tmp;
				else
					return x = tmp;
			*/

			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, double y) const {

			x = fmod (y, double(modulus));

			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, unsigned long y) const {

			x = fmod (float(y), modulus);

			if (x < 0) x += modulus;
			return x;
		}
		inline Element& init(Element& x, int y) const {

			x = fmod (float(y), modulus);

			if (x < 0) x += modulus;
			return x;
		}



		inline Element& assign(Element& x, const Element& y) const {
			return x = y;
		}


		inline bool areEqual (const Element &x, const Element &y) const {
			return x == y;
		}

		inline  bool isZero (const Element &x) const {
			return x == 0.;
		}

		inline bool isOne (const Element &x) const {
			return x == 1.;
		}

		inline Element &add (Element &x, const Element &y, const Element &z) const {
			x = y + z;
			if ( x >= modulus ) x -= modulus;
			return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const {
			x = y - z;
			if (x < 0) x += modulus;
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const {
			float tmp= y*z;
			x= fmod(tmp, modulus);
			//x= tmp - floor(tmp*inv_modulus)*modulus;

			return x;
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const {
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		inline Element &neg (Element &x, const Element &y) const {
			if(y == 0) return x = 0;
			else return x = modulus - y;
		}

		inline Element &inv (Element &x, const Element &y) const {
			// The extended Euclidean algoritm
			int x_int, y_int, q, tx, ty, temp;
			x_int = int (modulus);
			y_int = int (y);
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}

			if (tx < 0) tx += (int)modulus;

			// now x_int = gcd (modulus,residue)
			return x = (float)tx;


		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const {
			float tmp = a * x + y;
			return r= fmod(tmp, modulus);
			//return r= tmp- floor(tmp*inv_modulus)*modulus;

		}

		inline Element &addin (Element &x, const Element &y) const {
			x += y;
			if (  x >= modulus ) x -= modulus;
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const {
			x -= y;
			if (x < 0.) x += modulus;
			return x;
		}

		inline Element &mulin (Element &x, const Element &y) const {
			return mul(x,x,y);
		}

		inline Element &divin (Element &x, const Element &y) const {
			return div(x,x,y);
		}

		inline Element &negin (Element &x) const {
			if (x == 0.) return x;
			else return x = modulus - x;
		}

		inline Element &invin (Element &x) const {
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
			float tmp = r + a * x;
			return r = fmod(tmp, modulus);

			//return r= tmp- floor(tmp*inv_modulus)*modulus;
		}

		static inline float getMaxModulus()
			{ return 4096.0; } // floor( 2^12 )

	};

	template <>
	class FieldAXPY<Modular<float> > {
	public:

		typedef float Element;
		typedef Modular<float> Field;

		FieldAXPY (const Field &F) : _F (F) , //_invmod(1./_F.modulus),
					     _y(0.) , _bound( (float) ( (1UL << 23) - (int) (_F.modulus*_F.modulus))) {}

		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F),// _invmod(faxpy._invmod) ,
		_y(faxpy._y), _bound(faxpy._bound) {}

		FieldAXPY<Modular<float> > &operator = (const FieldAXPY &faxpy) {
			_F = faxpy._F;
			//_invmod= faxpy._invmod;
			_y= faxpy._y;
			_bound= faxpy._bound;
			return *this;
		}

            inline Element& mulacc (const Element &a, const Element &x) {
                Element tmp= a*x;
                return accumulate(tmp);
            }

            inline Element& accumulate (const Element &tmp) {
                _y += tmp;
                if (_y > _bound)
                    return _y = fmod (_y, _F.modulus);
                else
                    return _y;
            }

		inline Element& get (Element &y) {
			_y = fmod (_y, _F.modulus);
			return y=_y ;
		}

		inline FieldAXPY &assign (const Element y) {
			_y = y;
			return *this;
		}

		inline void reset() {
			_y = 0.;
		}

	private:

		Field _F;
		//float _invmod;
		float _y;
		float _bound;
	};


	template <>
	class DotProductDomain<Modular<float> > : private virtual VectorDomainBase<Modular<float> > {
	private:
		float _bound;
		size_t _nmax;
		//float _invmod;

	public:
		typedef float Element;
		DotProductDomain (const Modular<float> &F)
			: VectorDomainBase<Modular<float> > (F), _bound( (float) ( (1<<23) - (int) (_F.modulus*_F.modulus)))//, _invmod(1./_F.modulus)
			{
				_nmax= (size_t)floor((float(1<<11)* float(1<<12))/ (_F.modulus * _F.modulus));
			}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const {

			float y = 0.;
			float t = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y = fmod(y, _F.modulus);
			}
			else{
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmod(y, _F.modulus);
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmod(y, _F.modulus);
				y = fmod(t, _F.modulus);
			}
			return res = y;
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const {

			float y = 0.;
			float t =0.;


			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmod(y, _F.modulus);
			}
			else {
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmod(y, _F.modulus);
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmod(y, _F.modulus);
				y = fmod(t, _F.modulus);
			}
			return res = y;
		}
	};
}

#include "linbox/randiter/modular.h"

#endif
