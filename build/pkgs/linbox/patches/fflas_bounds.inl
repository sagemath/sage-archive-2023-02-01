/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* fflas/fflas_bounds.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifdef _LINBOX_LINBOX_CONFIG_H
#define FFLAS_INT_TYPE Integer
#else
#define FFLAS_INT_TYPE long unsigned int
#endif

/* From M Abshoff */
#if defined(__sun) && defined(__SunOS_5_9)
#define lround(x) my_lround(x)
static long my_lround(double x)
{
       return (long) ((x) >= 0 ? (x) + 0.5 : (x) - 0.5);
}
#endif


//---------------------------------------------------------------------
// DotProdBound :
// Computes the maximal dimension k so that the product A*B+beta.C over Z,
// where A is m*k and B is k*n can be performed correctly with w Winograd
// recursion levels on the number of bits of the floating point mantissa
//---------------------------------------------------------------------
template  <class Field>
inline size_t FFLAS::DotProdBoundCompute (const Field& F, const size_t w,
					  const typename Field::Element& beta,
					  const FFLAS_BASE base){

	typename Field::Element mone;
	static FFLAS_INT_TYPE p;
	F.characteristic(p);
	F.init (mone, -1.0);
	size_t kmax;

	unsigned long mantissa = (base == FflasDouble)? DOUBLE_MANTISSA : FLOAT_MANTISSA;

	if (p == 0)
		kmax = 2;
	else
		if (w > 0) {
			size_t ex=1;
			for (size_t i=0; i < w; ++i) 	ex *= 3;
			//FFLAS_INT_TYPE c = (p-1)*(ex)/2; //bound for a centered representation
			long long c;
#ifndef _LINBOX_LINBOX_CONFIG_H
			if (F.balanced)
				c = (p-1)*(ex)/2; // balanced representation
			else
#endif
				c = (p-1)*(1+ex)/2;
			kmax =  lround(( double(1ULL << mantissa) /double(c*c) + 1)*(1ULL << w));
			if (kmax ==  ( 1ULL << w))
				kmax = 2;
		}
		else{
			long long c = p-1;
			long long cplt=0;
			if (!F.isZero (beta))
				if (F.isOne (beta) || F.areEqual (beta, mone))
					cplt = c;
				else cplt = c*c;
			kmax =  lround(( double((1ULL << mantissa) - cplt)) /double(c*c));
			if (kmax  < 2)
				kmax = 2;
		}
	//cerr<<"kmax = "<<kmax<<endl;
	return  (size_t) MIN(kmax,1ULL<<31);
}



template  < class Field >
inline size_t FFLAS::DotProdBound (const Field& F, const size_t w,
				   const typename Field::Element& beta,
				   const FFLAS_BASE base)
{
	static Field G = F;
	static FFLAS_INT_TYPE pig;
	G.characteristic(pig);
	FFLAS_INT_TYPE pif;
	F.characteristic(pif);
	static typename Field::Element b = beta;
	static size_t w2 = w;
	static FFLAS_BASE base2 = base;
	static size_t kmax = DotProdBoundCompute (F, w, beta, base);
     	if ((b != beta) || (pif != pig) ||  (w2 != w) || (base2!=base)) {
		G = F;
		w2 = w;
		b = beta;
		base2 = base;
		kmax =  DotProdBoundCompute (F, w, beta, base);
	}
	return kmax;
}

//---------------------------------------------------------------------
// TRSMBound
// Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^53
//---------------------------------------------------------------------
size_t bound_compute_double(const long long pi) {

	long long p=pi,p1=1,p2=1;
	size_t nmax=1;
	double max = ( (  1ULL<<(DOUBLE_MANTISSA+1) )/(p-1));
	while ( (p1 + p2) < max ){
		p1*=p;
		p2*=p-2;
		nmax++;
	}
	nmax--;
	return nmax;
}
size_t bound_compute_double_balanced(const long long pi) {

	long long p=(pi+1)/2,p1=1;
	size_t nmax=0;
	double max = ( (  1ULL<<(DOUBLE_MANTISSA))/(p-1));
	while ( (p1) < max ){
		p1*=p;
		nmax++;
	}
	return nmax;
}
size_t bound_compute_float(const long long pi) {

	long long p=pi,p1=1,p2=1;
	size_t nmax=1;
	double max = ( (  1ULL<<(FLOAT_MANTISSA+1) )/(p-1));
	while ( (p1 + p2) < max ){
		p1*=p;
		p2*=p-2;
		nmax++;
	}
	nmax--;
	return nmax;
}
size_t bound_compute_float_balanced(const long long pi) {

	long long p=(pi+1)/2,p1=1;
	size_t nmax=0;
	double max = ( (  1ULL<<(FLOAT_MANTISSA))/(p-1));
	while ( (p1) < max ){
		p1*=p;
		nmax++;
	}
	return nmax;
}

template <class Field>
inline size_t
FFLAS::TRSMBound (const Field& F) {
	return callTRSMBound<typename Field::Element> () (F);
}

template<class Element>
class FFLAS::callTRSMBound {
public:
	template <class Field>
	size_t operator () (const Field& F) {
		FFLAS_INT_TYPE pi;
		F.characteristic(pi);
		static long unsigned int p=pi;
		static size_t nmax=bound_compute_double(p);
		if (p == pi)
			return nmax;
		else
			return nmax=bound_compute_double(p=pi);
	}
};

template<>
class FFLAS::callTRSMBound<double> {
public:
	template <class Field>
	size_t operator () (const Field& F) {
		FFLAS_INT_TYPE pi;
		F.characteristic(pi);
		static FFLAS_INT_TYPE p=pi;
#ifdef _LINBOX_LINBOX_CONFIG_H
		static size_t nmax = bound_compute_double(pi);
#else
		static size_t nmax = (F.balanced) ? bound_compute_double_balanced(pi) : bound_compute_double(pi);
#endif
		if (p == pi)
			return nmax;
		else
#ifdef _LINBOX_LINBOX_CONFIG_H
			return nmax= bound_compute_double (p=pi); //(F.balanced) ? bound_compute_balanced(p=pi) : bound_compute(p=pi);
#else
        	return (F.balanced) ? bound_compute_double_balanced(p=pi) : bound_compute_double(p=pi);
#endif
	}
};

template<>
class FFLAS::callTRSMBound<float> {
public:
	template <class Field>
	size_t operator () (const Field& F) {
		FFLAS_INT_TYPE pi;
		F.characteristic(pi);
		static FFLAS_INT_TYPE p=pi;
#ifdef _LINBOX_LINBOX_CONFIG_H
		static size_t nmax = bound_compute_float(pi);
#else
		static size_t nmax = (F.balanced) ? bound_compute_float_balanced(pi) : bound_compute_float(pi);
#endif
		if (p == pi)
			return nmax;
		else
#ifdef _LINBOX_LINBOX_CONFIG_H
			return nmax= bound_compute_float (p=pi); //(F.balanced) ? bound_compute_balanced(p=pi) : bound_compute(p=pi);
#else
        	return (F.balanced) ? bound_compute_float_balanced(p=pi) : bound_compute_float(p=pi);
#endif
	}
};



/** Automatic computation of the parameters for Matrix Multiplication
 *
 */

// To be tuned up: use statics to perform the computations only once for a given set of parameters
template <class Element>
class FFLAS::setMatMulParam {
 public:
	template <class Field>
		void operator() (const Field& F, const size_t m,
				 const typename Field::Element& beta,
				 size_t& w, FFLAS_BASE& base, size_t& kmax){

		// Strategy : determine Winograd's recursion first, then choose appropriate
		// floating point representation, and finally the blocking dimension.
		// Can be improved for some cases.
		static size_t ms = m;
		static size_t ks = WinoSteps (m);

		if (m != ms)
			ks = WinoSteps (ms = m);

		w = ks;

		base = BaseCompute<typename Field::Element> ()(F, w);
		kmax = DotProdBound (F, w, beta, base);
	}
};

template<>
class FFLAS::setMatMulParam<double> {
public:
	template <class Field>
	void operator () (const Field& F, const size_t m,
			  const typename Field::Element& beta,
			  size_t& w, FFLAS_BASE& base, size_t& kmax){

		static size_t ms = m;
		static size_t ks = WinoSteps (m);

		if (m != ms)
			ks = WinoSteps (ms = m);

		w = ks;
		base = FflasDouble;
		kmax = DotProdBound (F, w, beta, base);
	}
};

template<>
class FFLAS::setMatMulParam<float> {
public:
	template <class Field>
	void operator () (const Field& F, const size_t m,
			  const typename Field::Element& beta,
			  size_t& w, FFLAS_BASE& base, size_t& kmax){

		static size_t ms = m;
		static size_t ks = WinoSteps (m);

		if (m != ms)
			ks = WinoSteps (ms = m);

		w = ks;
		base = FflasFloat;
		kmax = DotProdBound (F, w, beta, base);
	}
};

inline size_t FFLAS::WinoSteps (const size_t m) {

	size_t w=0;
	size_t mt = m;
	while (mt >= WINOTHRESHOLD){
		w++;
		mt/=2;
	}
	//cerr<<"Winostep = "<<w<<endl;
	return w;
}


template <class Element>
class FFLAS::BaseCompute {
public:
	template <class Field>
	FFLAS::FFLAS_BASE operator() (const Field& F, const size_t w){

		FFLAS_INT_TYPE pi;
		F.characteristic(pi);
		FFLAS_BASE base;
		switch (w) {
		case 0: base = (pi < FLOAT_DOUBLE_THRESHOLD_0)? FflasFloat : FflasDouble;
			break;
		case 1:  base = (pi < FLOAT_DOUBLE_THRESHOLD_1)? FflasFloat : FflasDouble;
			break;
		case 2:  base = (pi < FLOAT_DOUBLE_THRESHOLD_2)? FflasFloat : FflasDouble;
			break;
		default: base = FflasDouble;
			break;
		}
		//cerr<<"base = "<<base<<endl;
		return base;
		// 	case 0: return (pi < FLOAT_DOUBLE_THRESHOLD_0)? FflasFloat : FflasDouble;
		// 	case 1: return (pi < FLOAT_DOUBLE_THRESHOLD_1)? FflasFloat : FflasDouble;
		// 	case 2: return (pi < FLOAT_DOUBLE_THRESHOLD_2)? FflasFloat : FflasDouble;
		// 	default: return FflasDouble;
	}
};

template<>
class FFLAS::BaseCompute<double> {
public:
	template <class Field>
	FFLAS::FFLAS_BASE operator() (const Field& F, const size_t w){
		return FflasDouble;
	}
};

template<>
class FFLAS::BaseCompute<float> {
public:
	template <class Field>
	FFLAS::FFLAS_BASE operator() (const Field& F, const size_t w){
		return FflasFloat;
	}
};
