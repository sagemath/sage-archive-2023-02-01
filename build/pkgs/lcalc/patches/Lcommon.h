// When MPFR is enabled and double is passed to a templated function
// The function should use precise(ttype) to make sure calculations run
// within the more precise type
#define precise(T) typename precise<T>::precise_type
template<class T> struct precise { typedef T precise_type; };
template<> struct precise<double> { typedef Double precise_type; };
template<> struct precise<complex<double> > { typedef Complex precise_type; };
#ifdef USE_MPFR
template<> struct precise<long long> { typedef double precise_type; };
#else
template<> struct precise<long long> { typedef long long precise_type; };
#endif
//typedef precise<long long>::precise_type Long;
typedef long long Long;


// To aid conversion to int value
inline double lcalc_to_double(const Double& x) {
#ifdef USE_MPFR
	return x.get_d();
#else
	return x;
#endif
}
#ifdef USE_MPFR
inline double lcalc_to_double(const double& x) { return x; }
#endif
inline double lcalc_to_double(const long double& x) { return x; }
inline double lcalc_to_double(const int& x) { return x; }
inline double lcalc_to_double(const long long& x) { return x; }
inline double lcalc_to_double(const short& x) { return x; }
inline double lcalc_to_double(const char& x) { return x; }
inline double lcalc_to_double(const long int& x) { return x; }
inline double lcalc_to_double(const unsigned int& x) { return x; }
inline double lcalc_to_double(const long unsigned int& x) { return x; }
#define Int(x) (int)(lcalc_to_double(x))
#define Long(x) (Long)(lcalc_to_double(x))
#define double(x) (double)(lcalc_to_double(x))

template<class T> inline int sn(T x)
{
 if (x>=0) return 1;
 else return -1;
}


const bool outputSeries=true;	  // Whether to output the coefficients or just the answer

// Loop i from m to n
// Useful in tidying up most for loops
#define loop(i,m,n) for(typeof(m) i=(m); i!=(n); i++)

// A class for calculations involving polynomials of small degree
// Not efficient enough for huge polynomials
//
// Polynomial of fixed degree N-1 with coefficients over T
template<class T=Complex> struct smallPoly {
	valarray<T> coeffs;
	typedef smallPoly<T> poly;
	int N;
	T zero;

	smallPoly(int N=0) {
		this->N=N;
		coeffs.resize(N);
		zero=0;
	}

	void resize(int N) {
		coeffs.resize(N);
		loop(i,this->N,N)
			coeffs[i]=zero;
		this->N=N;
	}

	// Access a coefficient as a subscript
	T& operator[](int n) {
		return (n<0 ? zero : coeffs[n]);
	}
	const T& operator[](int n) const {
		return (n<0 ? zero : coeffs[n]);
	}

	// Multiply two polys together, truncating the result
	friend poly operator*(const poly& f, const poly& g) {
		poly result(max(f.N,g.N));
		loop(i,0,f.N) result[i]=0;
		double test;
		loop(i,0,f.N) loop(j,0,result.N-i) {
			result[i+j]+=f[i]*g[j];
		}
		return result;
	}

	// Divide (in place) by (x-1/a) = (1-ax)/(-a) with a!=0
	void divideSpecial(T a) {
		loop(i,1,N) coeffs[i]+=coeffs[i-1]*a;
		coeffs*= -a;
	}

	// Evaluate the polynomial at x
	template<class U> U operator()(U x) {
		U sum=coeffs[0];
		U cur=1;

		loop(i,1,N) {
			cur*=x;
			sum+=coeffs[i]*x;
		}
		return sum;
	}

	// Output to stdout
	friend ostream& operator<<(ostream& s, const poly p) {
		loop(i,0,p.N) s << (i ? " + " : "") << p[i] << "x^" << i;
		return s;
	}

	// Arithmetic
	template<class U> poly& operator*=(const U& x) {
		coeffs*=x;
		return *this;
	}
	template<class U> poly& operator/=(const U& x) {
		coeffs/=x;
		return *this;
	}
};

