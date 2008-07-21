/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* ffpack.h
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __FFPACK_H
#define __FFPACK_H

#ifdef _LINBOX_LINBOX_CONFIG_H
#include "linbox/fflas/fflas.h"
#else
#include "fflas.h"
#endif

#include <list>
#include <vector>

#ifdef _LINBOX_LINBOX_CONFIG_H
namespace LinBox{
#endif

// The use of the small size LQUP is currently disabled:
// need for a better handling of element base (double, float, generic) combined
// with different thresholds.
// TransPosed version has to be implemented too.
#define __FFPACK_LUDIVINE_CUTOFF 0

#define __FFPACK_CHARPOLY_THRESHOLD 30

	/**
	 * \brief Set of elimination based routines for dense linear algebra
	 * with matrices over finite prime field of characteristic less than 2^26.
	 *
	 *  This class only provides a set of static member functions.
	 *  No instantiation is allowed.
	 *
	 * It enlarges the set of BLAS routines of the class FFLAS, with higher
	 * level routines based on elimination.
	 \ingroup ffpack
	 */

class FFPACK : public FFLAS {


public:
	enum FFPACK_LUDIVINE_TAG { FfpackLQUP=1,
				   FfpackSingular=2};

	enum FFPACK_CHARPOLY_TAG { FfpackLUK=1,
				   FfpackKG=2,
				   FfpackHybrid=3,
				   FfpackKGFast=4,
				   FfpackDanilevski=5,
				   FfpackArithProg=6,
				   FfpackKGFastG=7};

	class CharpolyFailed{};

	enum FFPACK_MINPOLY_TAG { FfpackDense=1,
				  FfpackKGF=2};

	/**
	 * Computes the rank of the given matrix using a LQUP factorization.
	 * The input matrix is modified.
	 * @param M row dimension of the matrix
	 * @param N column dimension of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 */
	/// using LQUP factorization.
	template <class Field>
	static size_t
	Rank( const Field& F, const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda)
	{
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		size_t R = LUdivine (F, FflasNonUnit, FflasNoTrans, M, N,
				     A, lda, P, Q, FfpackLQUP);
		delete[] Q;
		delete[] P;
		return R;
 	}

	/**
	 * Returns true if the given matrix is singular.
	 * The method is a block elimination with early termination
	 * The input matrix is modified.
	 * @param M row dimension of the matrix
	 * @param N column dimension of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 */
	/// using LQUP factorization  with early termination.
	template <class Field>
	static bool
	IsSingular( const Field& F, const size_t M, const size_t N,
		    typename Field::Element * A, const size_t lda)
	{
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		bool singular  = !LUdivine (F, FflasNonUnit, FflasNoTrans, M, N,
					    A, lda, P, Q, FfpackSingular);

		delete[] P;
		delete[] Q;
		return singular;
 	}

	/**
	 * Returns the determinant of the given matrix.
	 * The method is a block elimination with early termination
	 * The input matrix is modified.
	 * @param M row dimension of the matrix
	 * @param N column dimension of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 */
	///  using LQUP factorization  with early termination.
	template <class Field>
	static typename Field::Element
	Det( const Field& F, const size_t M, const size_t N,
	     typename Field::Element * A, const size_t lda){

		typename Field::Element det;
		bool singular;
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		singular  = !LUdivine (F, FflasNonUnit, FflasNoTrans,  M, N,
				       A, lda, P, Q, FfpackSingular);
		if (singular){
			F.init(det,0.0);
			delete[] P;
			delete[] Q;
			return det;
		}
		else{
			F.init(det,1.0);
			typename Field::Element *Ai=A;
			for (; Ai < A+ M*lda+N; Ai+=lda+1 )
				F.mulin( det, *Ai );
			int count=0;
			for (size_t i=0;i<N;++i)
				if (P[i] != i) ++count;

			if ((count&1) == 1)
				F.negin(det);
		}
		delete[] P;
		delete[] Q;
		return det;
 	}

	/**
	 * Solve the system A X = B or X A = B, using the LQUP decomposition of A
	 * already computed inplace with LUdivine(FflasNoTrans, FflasNonUnit).
	 * Version for A square.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 *
	 * @param Side Determine wheter the resolution is left or right looking.
	 * @param M row dimension of B
	 * @param N col dimension of B
	 * @param R rank of A
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param P column permutation of the LQUP decomposition of A
 	 * @param Q column permutation of the LQUP decomposition of A
	 * @param B Right/Left hand side matrix. Initially stores B, finally stores the solution X.
	 * @param ldb leading dimension of B
	 * @info Succes of the computation: 0 if successfull, >0 if system is inconsistent
	 */
	template <class Field>
	static void
	fgetrs (const Field& F,
		const FFLAS_SIDE Side,
		const size_t M, const size_t N, const size_t R,
		typename Field::Element *A, const size_t lda,
		const size_t *P, const size_t *Q,
		typename Field::Element *B, const size_t ldb,
		int * info){

		static typename Field::Element zero, one, mone;
		F.init (zero, 0.0);
		F.init (one, 1.0);
		F.neg(mone, one);
		*info =0;
		if (Side == FflasLeft) { // Left looking solve A X = B

			solveLB2 (F, FflasLeft, M, N, R, A, lda, Q, B, ldb);

			applyP (F, FflasLeft, FflasNoTrans, N, 0, R, B, ldb, Q);

			bool consistent = true;
			for (size_t i = R; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					if (!F.isZero (*(B + i*ldb + j)))
						consistent = false;
			if (!consistent) {
				std::cerr<<"System is inconsistent"<<std::endl;
				*info = 1;
			}
			// The last rows of B are now supposed to be 0
			//			for (size_t i = R; i < M; ++i)
			// 				for (size_t j = 0; j < N; ++j)
			// 					*(B + i*ldb + j) = zero;

			ftrsm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit,
			       R, N, one, A, lda , B, ldb);

			applyP (F, FflasLeft, FflasTrans, N, 0, R, B, ldb, P);

		} else { // Right Looking X A = B

			applyP (F, FflasRight, FflasTrans, M, 0, R, B, ldb, P);

			ftrsm (F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit,
			       M, R, one, A, lda , B, ldb);

			fgemm (F, FflasNoTrans, FflasNoTrans, M, N-R, R, one,
			       B, ldb, A+R, lda, mone, B+R, ldb);

			bool consistent = true;
			for (size_t i = 0; i < M; ++i)
				for (size_t j = R; j < N; ++j)
					if (!F.isZero (*(B + i*ldb + j)))
						consistent = false;
			if (!consistent) {
				std::cerr<<"System is inconsistent"<<std::endl;
				*info = 1;
			}
			// The last cols of B are now supposed to be 0

			applyP (F, FflasRight, FflasNoTrans, M, 0, R, B, ldb, Q);

			solveLB2 (F, FflasRight, M, N, R, A, lda, Q, B, ldb);
		}
	}

	/**
	 * Solve the system A X = B or X A = B, using the LQUP decomposition of A
	 * already computed inplace with LUdivine(FflasNoTrans, FflasNonUnit).
	 * Version for A rectangular.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 *
	 * @param Side Determine wheter the resolution is left or right looking.
	 * @param M row dimension of A
	 * @param N col dimension of A
	 * @param NRHS number of columns (if Side = FflasLeft) or row (if Side = FflasRight) of the matrices X and B
	 * @param R rank of A
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param P column permutation of the LQUP decomposition of A
 	 * @param Q column permutation of the LQUP decomposition of A
	 * @param X solution matrix
	 * @param ldx leading dimension of X
	 * @param B Right/Left hand side matrix.
	 * @param ldb leading dimension of B
	 * @param info Succes of the computation: 0 if successfull, >0 if system is inconsistent
	 */
	template <class Field>
	static typename Field::Element *
	fgetrs (const Field& F,
		const FFLAS_SIDE Side,
		const size_t M, const size_t N, const size_t NRHS, const size_t R,
		typename Field::Element *A, const size_t lda,
		const size_t *P, const size_t *Q,
		typename Field::Element *X, const size_t ldx,
		const typename Field::Element *B, const size_t ldb,
		int * info) {

		static typename Field::Element zero, one, mone;
		F.init (zero, 0.0);
		F.init (one, 1.0);
		F.neg(mone, one);
		*info =0;

		typename Field::Element* W;
		size_t ldw;

		if (Side == FflasLeft) { // Left looking solve A X = B

			// Initializing X to 0 (to be optimized)
			for (size_t i = 0; i <N; ++i)
				for (size_t j=0; j< NRHS; ++j)
					F.assign (*(X+i*ldx+j), zero);

			if (M > N){ // Cannot copy B into X
				W = new typename Field::Element [M*NRHS];
				ldw = NRHS;
				for (size_t i=0; i < M; ++i)
					fcopy (F, NRHS, W + i*ldw, 1, B + i*ldb, 1);

				solveLB2 (F, FflasLeft, M, NRHS, R, A, lda, Q, W, ldw);

				applyP (F, FflasLeft, FflasNoTrans, NRHS, 0, R, W, ldw, Q);

				bool consistent = true;
				for (size_t i = R; i < M; ++i)
					for (size_t j = 0; j < NRHS; ++j)
						if (!F.isZero (*(W + i*ldw + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					delete[] W;
					return X;
				}
				// Here the last rows of W are supposed to be 0

				ftrsm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit,
				       R, NRHS, one, A, lda , W, ldw);

				for (size_t i=0; i < R; ++i)
					fcopy (F, NRHS, X + i*ldx, 1, W + i*ldw, 1);

				delete[] W;
				applyP (F, FflasLeft, FflasTrans, NRHS, 0, R, X, ldx, P);

			} else { // Copy B to X directly
				for (size_t i=0; i < M; ++i)
					fcopy (F, NRHS, X + i*ldx, 1, B + i*ldb, 1);

				solveLB2 (F, FflasLeft, M, NRHS, R, A, lda, Q, X, ldx);

				applyP (F, FflasLeft, FflasNoTrans, NRHS, 0, R, X, ldx, Q);

				bool consistent = true;
				for (size_t i = R; i < M; ++i)
					for (size_t j = 0; j < NRHS; ++j)
						if (!F.isZero (*(X + i*ldx + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					return X;
				}
				// Here the last rows of W are supposed to be 0

				ftrsm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit,
				       R, NRHS, one, A, lda , X, ldx);

				applyP (F, FflasLeft, FflasTrans, NRHS, 0, R, X, ldx, P);
			}
			return X;

		} else { // Right Looking X A = B

			for (size_t i = 0; i <NRHS; ++i)
				for (size_t j=0; j< M; ++j)
					F.assign (*(X+i*ldx+j), zero);

			if (M < N) {
				W = new typename Field::Element [NRHS*N];
				ldw = N;
				for (size_t i=0; i < NRHS; ++i)
					fcopy (F, N, W + i*ldw, 1, B + i*ldb, 1);

				applyP (F, FflasRight, FflasTrans, NRHS, 0, R, W, ldw, P);

				ftrsm (F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit,
				       NRHS, R, one, A, lda , W, ldw);

				fgemm (F, FflasNoTrans, FflasNoTrans, NRHS, N-R, R, one,
				       W, ldw, A+R, lda, mone, W+R, ldw);

				bool consistent = true;
				for (size_t i = 0; i < NRHS; ++i)
					for (size_t j = R; j < N; ++j)
						if (!F.isZero (*(W + i*ldw + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					delete[] W;
					return X;
				}
				// The last N-R cols of W are now supposed to be 0
				for (size_t i=0; i < NRHS; ++i)
					fcopy (F, R, X + i*ldx, 1, W + i*ldb, 1);
				delete[] W;
				applyP (F, FflasRight, FflasNoTrans, NRHS, 0, R, X, ldx, Q);

				solveLB2 (F, FflasRight, NRHS, M, R, A, lda, Q, X, ldx);

			} else {
				for (size_t i=0; i < NRHS; ++i)
					fcopy (F, N, X + i*ldx, 1, B + i*ldb, 1);

				applyP (F, FflasRight, FflasTrans, NRHS, 0, R, X, ldx, P);

				ftrsm (F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit,
				       NRHS, R, one, A, lda , X, ldx);

				fgemm (F, FflasNoTrans, FflasNoTrans, NRHS, N-R, R, one,
				       X, ldx, A+R, lda, mone, X+R, ldx);

				bool consistent = true;
				for (size_t i = 0; i < NRHS; ++i)
					for (size_t j = R; j < N; ++j)
						if (!F.isZero (*(X + i*ldx + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					return X;
				}
				// The last N-R cols of W are now supposed to be 0

				applyP (F, FflasRight, FflasNoTrans, NRHS, 0, R, X, ldx, Q);

				solveLB2 (F, FflasRight, NRHS, M, R, A, lda, Q, X, ldx);

			}
			return X;
		}
	}
	/**
	 * @brief Square system solver
	 * @param Field The computation domain
	 * @param Side Determine wheter the resolution is left or right looking
	 * @param M row dimension of B
	 * @param N col dimension of B
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param P column permutation of the LQUP decomposition of A
 	 * @param Q column permutation of the LQUP decomposition of A
	 * @param B Right/Left hand side matrix. Initially contains B, finally contains the solution X.
	 * @param ldb leading dimension of B
	 * @param info Success of the computation: 0 if successfull, >0 if system is inconsistent
 	 * @return the rank of the system
 	 *
	 * Solve the system A X = B or X A = B.
	 * Version for A square.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 */
	template <class Field>
	static size_t
	fgesv (const Field& F,
	       const FFLAS_SIDE Side,
	       const size_t M, const size_t N,
	       typename Field::Element *A, const size_t lda,
	       typename Field::Element *B, const size_t ldb,
	       int * info){

		size_t Na;
		if (Side == FflasLeft)
			Na = M;
		else
			Na = N;

		size_t* P = new size_t[Na];
		size_t* Q = new size_t[Na];

		size_t R = LUdivine (F, FflasNonUnit, FflasNoTrans, Na, Na, A, lda, P, Q, FfpackLQUP);

		fgetrs (F, Side, M, N, R, A, lda, P, Q, B, ldb, info);

		delete[] P;
		delete[] Q;

		return R;
	}

	/**
	 * @brief Rectangular system solver
	 * @param Field The computation domain
	 * @param Side Determine wheter the resolution is left or right looking
	 * @param M row dimension of A
	 * @param N col dimension of A
	 * @param NRHS number of columns (if Side = FflasLeft) or row (if Side = FflasRight) of the matrices X and B
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param P column permutation of the LQUP decomposition of A
 	 * @param Q column permutation of the LQUP decomposition of A
	 * @param B Right/Left hand side matrix. Initially contains B, finally contains the solution X.
	 * @param ldb leading dimension of B
	 * @info Success of the computation: 0 if successfull, >0 if system is inconsistent
 	 * @return the rank of the system
 	 *
	 * Solve the system A X = B or X A = B.
	 * Version for A square.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 */
	template <class Field>
	static size_t
	fgesv (const Field& F,
	       const FFLAS_SIDE Side,
	       const size_t M, const size_t N, const size_t NRHS,
	       typename Field::Element *A, const size_t lda,
	       typename Field::Element *X, const size_t ldx,
	       const typename Field::Element *B, const size_t ldb,
	       int * info){

		size_t Nb,Mb;
		if (Side == FflasLeft){Nb = NRHS; Mb = N;}
		else {Nb = M; Mb = NRHS;}

		size_t* P = new size_t[N];
		size_t* Q = new size_t[M];

		size_t R = LUdivine (F, FflasNonUnit, FflasNoTrans, M, N, A, lda, P, Q, FfpackLQUP);

		fgetrs (F, Side, M, N, NRHS, R, A, lda, P, Q, X, ldx, B, ldb, info);

		delete[] P;
		delete[] Q;

		return R;
	}

	/**
	 * Solve the system Ax=b, using LQUP factorization and
	 * two triangular system resolutions.
	 * The input matrix is modified.
	 * @param M row dimension of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param x solution vector
	 * @param incX increment of x
	 * @param b right hand side vector
	 * @param incB increment of b
	 */
	/// Solve linear system using LQUP factorization.
	template <class Field>
	static typename Field::Element*
	Solve( const Field& F, const size_t M,
	       typename Field::Element * A, const size_t lda,
	       typename Field::Element * x, const int incx,
	       const typename Field::Element * b, const int incb ){
		typename Field::Element one, zero;
		F.init(one,1.0);
		F.init(zero,0.0);

		size_t *P = new size_t[M];
		size_t *rowP = new size_t[M];

		if (LUdivine( F, FflasNonUnit, FflasNoTrans, M, M, A, lda, P, rowP, FfpackLQUP) < M){
			std::cerr<<"SINGULAR MATRIX"<<std::endl;
			delete[] P;
			delete[] rowP;
			return x;
		}
		else{
			fcopy( F, M, x, incx, b, incb );

			ftrsv(F,  FflasLower, FflasNoTrans, FflasUnit, M,
			      A, lda , x, incx);
			ftrsv(F,  FflasUpper, FflasNoTrans, FflasNonUnit, M,
			      A, lda , x, incx);
			applyP( F, FflasRight, FflasTrans, M, 0, M, x, incx, P );
			delete[] rowP;
			delete[] P;

			return x;

		}
 	}

	/**
	 * Invert the given matrix in place
	 * or computes its nullity if it is singular.
	 * An inplace 2n^3 algorithm is used.
	 * @param M order of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param nullity dimension of the kernel of A
	 */
	/// Invert a matrix or return its nullity
	template <class Field>
	static typename Field::Element*
	Invert (const Field& F, const size_t M,
		typename Field::Element * A, const size_t lda,
		int& nullity){

		size_t * P = new size_t[M];
		size_t * Q = new size_t[M];
		size_t R =  ReducedColumnEchelonForm (F, M, M, A, lda, P, Q);
		nullity = M - R;
		applyP (F, FflasLeft, FflasTrans, M, 0, R, A, lda, P);
		delete [] P;
		delete [] Q;
		return A;
	}

	template <class Field>
	static typename Field::Element*
	Invert (const Field& F, const size_t M,
		typename Field::Element * A, const size_t lda,
		typename Field::Element * X, const size_t ldx,
		int& nullity){

		Invert (F,  M, A, lda, nullity);
		for (size_t i=0; i<M; ++i)
			fcopy (F, M, X+i*ldx, 1, A+i*lda,1);
		return X;

	}
	/**
	 * Invert the given matrix or computes its nullity if it is singular.
	 * An 2n^3 algorithm is used.
	 * The input matrix is modified.
	 * @param M order of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param X inverse of A
	 * @param ldx leading dimension of X
	 * @param nullity dimension of the kernel of A
	 */
	/// Invert a matrix or return its nullity
	template <class Field>
	static typename Field::Element*
	Invert2( const Field& F, const size_t M,
		 typename Field::Element * A, const size_t lda,
		 typename Field::Element * X, const size_t ldx,
		 int& nullity){

		typename Field::Element one, zero;
		F.init(one,1.0);
		F.init(zero,0.0);

		size_t *P = new size_t[M];
		size_t *rowP = new size_t[M];

		// Timer t1;
// 		t1.clear();
// 		t1.start();

		nullity = M - LUdivine( F, FflasNonUnit, FflasNoTrans, M, M, A, lda, P, rowP, FfpackLQUP);

// 		t1.stop();
		//cerr<<"LU --> "<<t1.usertime()<<endl;

		if (nullity > 0){
			delete[] P;
			delete[] rowP;
			return NULL;
		} else {
			// Initializing X to 0
// 			t1.clear();
// 			t1.start();
			for (size_t i=0; i<M; ++i)
				for (size_t j=0; j<M;++j)
					F.assign(*(X+i*ldx+j), zero);

			// X = L^-1 in n^3/3
 			ftrtri (F, FflasLower, FflasUnit, M, A, lda);
			for (size_t i=0; i<M; ++i){
				for (size_t j=i; j<M; ++j)
					F.assign(*(X +i*ldx+j),zero);
				F.assign (*(X+i*(ldx+1)), one);
			}
			for (size_t i=1; i<M; ++i)
 				fcopy (F, i, (X+i*ldx), 1, (A+i*lda), 1);
// 			t1.stop();
			//cerr<<"U^-1 --> "<<t1.usertime()<<endl;

			//invL( F, M, A, lda, X, ldx );
		       // X = Q^-1.X is not necessary since Q = Id

 			// X = U^-1.X
// 			t1.clear();
//			t1.start();
			ftrsm( F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit,
			       M, M, one, A, lda , X, ldx);
//			t1.stop();
			//cerr<<"ftrsm --> "<<t1.usertime()<<endl;

			// X = P^-1.X
			applyP( F, FflasLeft, FflasTrans, M, 0, M, X, ldx, P );

			delete[] P;
			delete[] rowP;
			return X;
		}
 	}

	/** RowRankProfile
	 * Computes the row rank profile of A.
	 *
	 * @param A: input matrix of dimension
	 * @param rklprofile: return the rank profile as an array of row indexes, of dimension r=rank(A)
	 *
	 * rkprofile is allocated during the computation.
	 * Returns R
	 */
	template <class Field>
	static size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
				      typename Field::Element* A, const size_t lda, size_t* rkprofile){
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		size_t R;

		R = LUdivine (F, FflasNonUnit, FflasNoTrans, M, N, A, lda, P, Q);
		rkprofile = new size_t[R];

		for (size_t i=0; i<R; ++i)
			rkprofile[i] = Q[i];
		delete[] P;
		delete[] Q;
		return R;
	}

	/** ColumnRankProfile
	 * Computes the column rank profile of A.
	 *
	 * @param A: input matrix of dimension
	 * @param rklprofile: return the rank profile as an array of row indexes, of dimension r=rank(A)
	 *
	 * A is modified
	 * rkprofile is allocated during the computation.
	 * Returns R
	 */
	template <class Field>
	static size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
					 typename Field::Element* A, const size_t lda, size_t* rkprofile){
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		size_t R;

		R = LUdivine (F, FflasNonUnit, FflasTrans, M, N, A, lda, P, Q);
		rkprofile = new size_t[R];

		for (size_t i=0; i<R; ++i)
			rkprofile[i] = Q[i];
		delete[] P;
		delete[] Q;
		return R;
	}

	/** RowRankProfileSubmatrixIndices
	 * Computes the indices of the submatrix r*r X of A whose rows correspond to
	 * the row rank profile of A.
	 *
	 * @param A: input matrix of dimension
	 * @param rowindices: array of the row indices of X in A
	 * @param colindices: array of the col indices of X in A
	 *
	 * rowindices and colindices are allocated during the computation.
	 * A is modified
	 * Returns R
	 */
	template <class Field>
	static size_t RowRankProfileSubmatrixIndices (const Field& F,
						      const size_t M, const size_t N,
						      typename Field::Element* A,
						      const size_t lda,
						      size_t*& rowindices,
						      size_t*& colindices,
						      size_t& R){
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];

		R = LUdivine (F, FflasNonUnit, FflasNoTrans, M, N, A, lda, P, Q);
		rowindices = new size_t[M];
		colindices = new size_t[N];
		for (size_t i=0; i<R; ++i){
			rowindices [i] = Q [i];
		}
		for (size_t i=0; i<N; ++i)
			colindices [i] = i;
		size_t tmp;
		for (size_t i=0; i<R; ++i){
			if (i != P[i]){
				tmp = colindices[i];
				colindices[i] = colindices[P[i]];
				colindices[P[i]] = tmp;
			}
		}

		delete[] P;
		delete[] Q;

		return R;
	}

	/** ColRankProfileSubmatrixIndices
	 * Computes the indices of the submatrix r*r X of A whose columns correspond to
	 * the column rank profile of A.
	 *
	 * @param A: input matrix of dimension
	 * @param rowindices: array of the row indices of X in A
	 * @param colindices: array of the col indices of X in A
	 *
	 * rowindices and colindices are allocated during the computation.
	 * A is modified
	 * Returns R
	 */
	template <class Field>
	static size_t ColRankProfileSubmatrixIndices (const Field& F,
						      const size_t M, const size_t N,
						      typename Field::Element* A,
						      const size_t lda,
						      size_t*& rowindices,
						      size_t*& colindices,
						      size_t& R){
		size_t *P = new size_t[M];
		size_t *Q = new size_t[N];

		R = LUdivine (F, FflasNonUnit, FflasTrans, M, N, A, lda, P, Q);
		rowindices = new size_t[M];
		colindices = new size_t[N];
		for (size_t i=0; i<R; ++i)
			colindices [i] = Q [i];

		for (size_t i=0; i<N; ++i)
			rowindices [i] = i;

		size_t tmp;
		for (size_t i=0; i<R; ++i){
			if (i != P[i]){
				tmp = rowindices[i];
				rowindices[i] = rowindices[P[i]];
				rowindices[P[i]] = tmp;
			}
		}
		delete[] P;
		delete[] Q;

		return R;
	}

	/** RowRankProfileSubmatrix
	 * Compute the r*r submatrix X of A, by picking the row rank profile rows of A
	 *
	 * @param A: input matrix of dimension M x N
	 * @param X: the output matrix
	 *
	 * A is not modified
	 * X is allocated during the computation.
	 * Returns R
	 */
	template <class Field>
	static size_t RowRankProfileSubmatrix (const Field& F,
					       const size_t M, const size_t N,
					       typename Field::Element* A,
					       const size_t lda,
					       typename Field::Element*& X, size_t& R){

		size_t * rowindices, * colindices;

		typename Field::Element * A2 = MatCopy (F, M, N, A, lda);

		RowRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

		X = new typename Field::Element[R*R];
		for (size_t i=0; i<R; ++i)
			for (size_t j=0; j<R; ++j)
				F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
		delete[] A2;
		delete[] rowindices;
		delete[] colindices;
		return R;
	}


	/** ColRankProfileSubmatrix
	 * Compute the r*r submatrix X of A, by picking the row rank profile rows of A
	 *
	 * @param A: input matrix of dimension M x N
	 * @param X: the output matrix
	 *
	 * A is not modified
	 * X is allocated during the computation.
	 * Returns R
	 */

	template <class Field>
	static size_t ColRankProfileSubmatrix (const Field& F, const size_t M, const size_t N,
					       typename Field::Element* A, const size_t lda,
					       typename Field::Element*& X, size_t& R){

		size_t * rowindices, * colindices;

		typename Field::Element * A2 = MatCopy (F, M, N, A, lda);

		ColRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

		X = new typename Field::Element[R*R];
		for (size_t i=0; i<R; ++i)
			for (size_t j=0; j<R; ++j)
				F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
		delete[] A2;
		delete[] colindices;
		delete[] rowindices;
		return R;
	}


	/**
	 * LQUPtoInverseOfFullRankMinor
	 * Suppose A has been factorized as L.Q.U.P, with rank r.
	 * Then Qt.A.Pt has an invertible leading principal r x r submatrix
	 * This procedure efficiently computes the inverse of this minor and puts it into X.
	 * NOTE: It changes the lower entries of A_factors in the process (NB: unless A was nonsingular and square)
	 *
	 * @param rank:       rank of the matrix.
	 * @param A_factors:  matrix containing the L and U entries of the factorization
	 * @param QtPointer:  theLQUP->getQ()->getPointer() (note: getQ returns Qt!)
	 * @param X:          desired location for output
	 */
	template <class Field>
	static typename Field::Element*
	LQUPtoInverseOfFullRankMinor( const Field& F, const size_t rank,
				      typename Field::Element * A_factors, const size_t lda,
				      const size_t* QtPointer,
				      typename Field::Element * X, const size_t ldx){

		typename Field::Element one, zero;
		F.init(one,1.0);
		F.init(zero,0.0);

		// upper entries are okay, just need to move up bottom ones
		const size_t* srcRow = QtPointer;
		for (size_t row=0; row<rank; row++, srcRow++)
			if (*srcRow != row) {
				typename Field::Element* oldRow = A_factors + (*srcRow) * lda;
				typename Field::Element* newRow = A_factors + row * lda;
				for (size_t col=0; col<row; col++, oldRow++, newRow++)
					F.assign(*newRow, *oldRow);
			}

		// X <- (Qt.L.Q)^(-1)
		//invL( F, rank, A_factors, lda, X, ldx);
		ftrtri (F, FflasLower, FflasUnit, rank, A_factors, lda);
		for (size_t i=0; i<rank; ++i)
			fcopy (F, rank, A_factors+i*lda, 1, X+i*ldx,1);

		// X = U^-1.X
		ftrsm( F, FflasLeft, FflasUpper, FflasNoTrans,
		       FflasNonUnit, rank, rank, one, A_factors, lda, X, ldx);

		return X;

 	}


	//---------------------------------------------------------------------
	// TURBO: rank computation algorithm
	//---------------------------------------------------------------------
	template <class Field>
	static size_t
	TURBO (const Field& F, const size_t M, const size_t N,
	       typename Field::Element* A, const size_t lda, size_t * P, size_t * Q, const size_t cutoff);

	/**
	 * Compute the LQUP factorization of the given matrix using
	 * a block agorithm and return its rank.
	 * The permutations P and Q are represented
	 * using LAPACK's convention.
	 * @param Diag  precise whether U should have a unit diagonal or not
	 * @param M matrix row dimension
	 * @param N matrix column dimension
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param P the column permutation
	 * @param Qt the transpose of the row permutation Q
	 * @param LuTag flag for setting the earling termination if the matrix
	 * is singular
	 */
	/// LQUP factorization.
	template <class Field>
	static size_t
	LUdivine (const Field& F, const FFLAS_DIAG Diag,  const FFLAS_TRANSPOSE trans,
		  const size_t M, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  size_t* P, size_t* Qt,
		  const FFPACK_LUDIVINE_TAG LuTag=FfpackLQUP,
		  const size_t cutoff=__FFPACK_LUDIVINE_CUTOFF);


	template<class Element>
	class callLUdivine_small;

	template <class Field>
	static size_t
	LUdivine_small (const Field& F, const FFLAS_DIAG Diag,  const FFLAS_TRANSPOSE trans,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			size_t* P, size_t* Q, const FFPACK_LUDIVINE_TAG LuTag=FfpackLQUP);

	template <class Field>
	static size_t
	LUdivine_gauss (const Field& F, const FFLAS_DIAG Diag,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			size_t* P, size_t* Q, const FFPACK_LUDIVINE_TAG LuTag=FfpackLQUP);



	/**
	 * Compute the inverse of a triangular matrix.
	 * @param Uplo whether the matrix is upper of lower triangular
	 * @param Diag whether the matrix if unit diagonal
	 *
	 */
	 template<class Field>
	 static void
	 ftrtri (const Field& F, const FFLAS_UPLO Uplo, const FFLAS_DIAG Diag,
		 const size_t N, typename Field::Element * A, const size_t lda){

		 static typename Field::Element one;
		 static typename Field::Element mone;
		 F.init(one,1.0);
		 F.init(mone,-1.0);
		 if (N == 1){
			 if (Diag == FflasNonUnit)
				 F.invin (*A);
		 } else {
			 size_t N1 = N/2;
			 size_t N2 = N - N1;
			 ftrtri (F, Uplo, Diag, N1, A, lda);
			 ftrtri (F, Uplo, Diag, N2, A + N1*(lda+1), lda);
			 if (Uplo == FflasUpper){
				 ftrmm (F, FflasLeft, Uplo, FflasNoTrans, Diag, N1, N2,
					one, A, lda, A + N1, lda);
				 ftrmm (F, FflasRight, Uplo, FflasNoTrans, Diag, N1, N2,
					mone, A + N1*(lda+1), lda, A + N1, lda);
			 } else {
				 ftrmm (F, FflasLeft, Uplo, FflasNoTrans, Diag, N2, N1,
					one, A + N1*(lda+1), lda, A + N1*lda, lda);
				 ftrmm (F, FflasRight, Uplo, FflasNoTrans, Diag, N2, N1,
					mone, A, lda, A + N1*lda, lda);
			 }
		 }
	 }


	/**
	 * Compute the product UL of the upper, resp lower triangular matrices U and L
	 * stored one above the other in the square matrix A.
	 * Diag == Unit if the matrix U is unit diagonal
	 *
	 */
	template<class Field>
	static void
	ftrtrm (const Field& F, const FFLAS_DIAG diag, const size_t N,
		typename Field::Element * A, const size_t lda){

		typename Field::Element one;
		F.init(one,1.0);

		if (N == 1)
			return;
		size_t N1 = N/2;
		size_t N2 = N-N1;

		ftrtrm (F, diag, N1, A, lda);

		fgemm (F, FflasNoTrans, FflasNoTrans, N1, N1, N2, one,
		       A+N1, lda, A+N1*lda, lda, one, A, lda);

		ftrmm (F, FflasRight, FflasLower, FflasNoTrans, (diag == FflasUnit) ? FflasNonUnit : FflasUnit, N1, N2, one, A + N1*(lda+1), lda, A + N1, lda);

		ftrmm (F, FflasLeft, FflasUpper, FflasNoTrans, diag, N2, N1, one, A + N1*(lda+1), lda, A + N1*lda, lda);

		ftrtrm (F, diag, N2, A + N1*(lda+1), lda);

	}

	/**
	 * Compute the Column Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ M \ V ] such that AU = C is a column echelon
	 * decomposition of A, with U = P^T [ V + Ir ] and C = M //+ Q [ Ir   ]
	 *                                  [ 0 In-r ]           //    [    0 ]
	 * Qt = Q^T
	 */
	template <class Field>
	static size_t
	ColumnEchelonForm (const Field& F, const size_t M, const size_t N,
			   typename Field::Element * A, const size_t lda,
			   size_t* P, size_t* Qt){

		typename Field::Element one, mone;
		F.init (one, 1.0);
		F.neg (mone, one);
		size_t r;

		// Timer t1;
// 		t1.clear();
// 		t1.start();
		r = LUdivine (F, FflasUnit, FflasNoTrans, M, N, A, lda, P, Qt);
		// t1.stop();
		//cerr<<"LU --> "<<t1.usertime()<<endl;

		// Timer t2;
// 		t2.clear();
// 		t2.start();
		ftrtri (F, FflasUpper, FflasUnit, r, A, lda);


		ftrmm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasUnit, r, N-r,
		       mone, A, lda, A+r, lda);

		// t2.stop();
		//cerr<<"U^-1 --> "<<t2.usertime()<<endl;

		return r;
	}

	/**
	 * Compute the Row Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ L \ M ] such that L A = R is a column echelon
	 * decomposition of A, with L =  [ L+Ir  0   ] P  and R = M
	 *                               [      In-r ]
	 * Qt = Q^T
	 */
	template <class Field>
	static size_t
	RowEchelonForm (const Field& F, const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			size_t* P, size_t* Qt){

		typename Field::Element one, mone;
		F.init (one, 1.0);
		F.neg (mone, one);
		size_t r;

		// Timer t1;
// 		t1.clear();
// 		t1.start();
		r = LUdivine (F, FflasUnit, FflasTrans,  M, N, A, lda, P, Qt);
		// t1.stop();
		//cerr<<"LU --> "<<t1.usertime()<<endl;

		// Timer t2;
// 		t2.clear();
// 		t2.start();
		ftrtri (F, FflasLower, FflasUnit, r, A, lda);


		ftrmm (F, FflasRight, FflasLower, FflasNoTrans, FflasUnit, M-r, r,
		       mone, A, lda, A+r*lda, lda);

		// t2.stop();
		//cerr<<"U^-1 --> "<<t2.usertime()<<endl;

		return r;
	}

	/**
	 * Compute the Reduced Column Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [  V  ] such that AU = R is a reduced column echelon
	 *                           [ M 0 ]
	 * decomposition of A, where U = P^T [ V      ] and R = Q [ Ir   ]
	 *                                   [ 0 In-r ]           [ M  0 ]
	 * Qt = Q^T
	 */
	template <class Field>
	static size_t
	ReducedColumnEchelonForm (const Field& F, const size_t M, const size_t N,
				  typename Field::Element * A, const size_t lda,
				  size_t* P, size_t* Qt){

		typename Field::Element one, mone;
		F.init (one, 1.0);
		F.neg (mone, one);
		size_t r;

		r = ColumnEchelonForm (F, M, N, A, lda, P, Qt);

// 		Timer t1;
// 		t1.clear();
// 		t1.start();

		// M = Q^T M
		for (size_t i=0; i<r; ++i){
			if ( Qt[i]> (size_t) i ){
				fswap( F, i+1,
				       A + Qt[i]*lda, 1,
				       A + i*lda, 1 );
			}
		}

		ftrtri (F, FflasLower, FflasNonUnit, r, A, lda);

		ftrmm (F, FflasRight, FflasLower, FflasNoTrans, FflasNonUnit, M-r, r,
		       one, A, lda, A+r*lda, lda);

		ftrtrm (F, FflasUnit, r, A, lda);
// 		t1.stop();
		//cerr<<"U^-1L^-1 --> "<<t1.usertime()<<endl;

		return r;

	}

	/**
	 * Compute the Reduced Row Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [  V  ] such that L A = R is a reduced row echelon
	 *                           [ M 0 ]
	 * decomposition of A, where L =  [ V      ] P^T and R =  [ Ir M  ] Q
	 *                                [ 0 In-r ]              [ 0     ]
	 * Qt = Q^T
	 */
	template <class Field>
	static size_t
	ReducedRowEchelonForm (const Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       size_t* P, size_t* Qt){

		typename Field::Element one, mone;
		F.init (one, 1.0);
		F.neg (mone, one);
		size_t r;

		r = RowEchelonForm (F, M, N, A, lda, P, Qt);

// 		Timer t1;
// 		t1.clear();
// 		t1.start();
		// M = M Q^T
		for (int i=0; i<r; ++i){
			if ( Qt[i]> (size_t) i ){
				fswap( F, i+1,
				       A + Qt[i], lda,
				       A + i, lda );
			}
		}

		ftrtri (F, FflasUpper, FflasNonUnit, r, A, lda);

		ftrmm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, r, N-r,
		       one, A, lda, A+r, lda);

		ftrtrm (F, FflasNonUnit, r, A, lda);

// 		t1.stop();
		//cerr<<"U^-1L^-1 --> "<<t1.usertime()<<endl;

		return r;

	}

	/** Apply a permutation submatrix of P (between ibeg and iend) to a matrix
	 * to (iend-ibeg) vectors of size M stored in A (as column for NoTrans
	 * and rows for Trans)
	 * Side==FflasLeft for row permutation Side==FflasRight for a column
	 * permutation
	 * Trans==FflasTrans for the inverse permutation of P
	 */
	template<class Field>
	static void
	applyP( const Field& F,
		const FFLAS_SIDE Side,
		const FFLAS_TRANSPOSE Trans,
		const size_t M, const int ibeg, const int iend,
		typename Field::Element * A, const size_t lda, const size_t * P ){

		if ( Side == FflasRight )
			if ( Trans == FflasTrans ){
				for ( size_t i=ibeg; i<(size_t) iend; ++i){
					if ( P[i]> i )
						fswap( F, M,
						       A + P[i]*1, lda,
						       A + i*1, lda );
				}
			}
			else{ // Trans == FflasNoTrans
				for (int i=iend-1; i>=ibeg; --i){
					if ( P[i]>(size_t)i ){
						fswap( F, M,
						       A + P[i]*1, lda,
						       A + i*1, lda );
					}
				}
			}
		else // Side == FflasLeft
			if ( Trans == FflasNoTrans ){
				for (size_t i=ibeg; i<(size_t)iend; ++i){
					if ( P[i]> (size_t) i )
						fswap( F, M,
						       A + P[i]*lda, 1,
						       A + i*lda, 1 );
				}
			}
			else{ // Trans == FflasTrans
				for (int i=iend-1; i>=ibeg; --i){
					if ( P[i]> (size_t) i ){
						fswap( F, M,
						       A + P[i]*lda, 1,
						       A + i*lda, 1 );
					}
				}
			}

	}

	/**
	 * Compute the characteristic polynomial of A using Krylov
	 * Method, and LUP factorization of the Krylov matrix
	 */
	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	CharPoly( const Field& F, std::list<Polynomial>& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  const FFPACK_CHARPOLY_TAG CharpTag= FfpackLUK);

	/**
	 * Compute the minimal polynomial of (A,v) using an LUP
	 * factorization of the Krylov Base (v, Av, .., A^kv)
	 * U,X must be (n+1)*n
	 * U contains the Krylov matrix and X, its LSP factorization
	 */
	template <class Field, class Polynomial>
	static Polynomial&
	MinPoly( const Field& F, Polynomial& minP, const size_t N,
		 const typename Field::Element *A, const size_t lda,
		 typename Field::Element* X, const size_t ldx, size_t* P,
		 const FFPACK_MINPOLY_TAG MinTag,
		 const size_t kg_mc, const size_t kg_mb, const size_t kg_j );


	// Solve L X = B or X L = B in place
	// L is M*M if Side == FflasLeft and N*N if Side == FflasRight, B is M*N.
	// Only the R non trivial column of L are stored in the M*R matrix L
	// Requirement :  so that L could  be expanded in-place
	template<class Field>
	static void
	solveLB( const Field& F, const FFLAS_SIDE Side,
		 const size_t M, const size_t N, const size_t R,
		 typename Field::Element * L, const size_t ldl,
		 const size_t * Q,
		 typename Field::Element * B, const size_t ldb ){

		typename Field::Element one, zero;
		F.init(one, 1.0);
		F.init(zero, 0.0);
		size_t LM = (Side == FflasRight)?N:M;
		for (int i=R-1; i>=0; --i){
			if (  Q[i] > (size_t) i){
				//for (size_t j=0; j<=Q[i]; ++j)
				//F.init( *(L+Q[i]+j*ldl), 0 );
				//std::cerr<<"1 deplacement "<<i<<"<-->"<<Q[i]<<endl;
				fcopy( F, LM-Q[i]-1, L+Q[i]*(ldl+1)+ldl,ldl, L+(Q[i]+1)*ldl+i, ldl );
				for ( size_t j=Q[i]*ldl; j<LM*ldl; j+=ldl)
					F.assign( *(L+i+j), zero );
			}
		}
		ftrsm( F, Side, FflasLower, FflasNoTrans, FflasUnit, M, N, one, L, ldl , B, ldb);
		//write_field(F,std::cerr<<"dans solveLB "<<endl,L,N,N,ldl);
		// Undo the permutation of L
		for (size_t i=0; i<R; ++i){
			if ( Q[i] > (size_t) i){
				//for (size_t j=0; j<=Q[i]; ++j)
				//F.init( *(L+Q[i]+j*ldl), 0 );
				fcopy( F, LM-Q[i]-1, L+(Q[i]+1)*ldl+i, ldl, L+Q[i]*(ldl+1)+ldl,ldl );
				for ( size_t j=Q[i]*ldl; j<LM*ldl; j+=ldl)
					F.assign( *(L+Q[i]+j), zero );
			}
		}
	}

	// Solve L X = B in place
	// L is M*M or N*N, B is M*N.
	// Only the R non trivial column of L are stored in the M*R matrix L
	template<class Field>
	static void
	solveLB2( const Field& F, const FFLAS_SIDE Side,
		  const size_t M, const size_t N, const size_t R,
		  typename Field::Element * L, const size_t ldl,
		  const size_t * Q,
		  typename Field::Element * B, const size_t ldb ){



		typename Field::Element Mone, one;
		F.init( Mone, -1.0 );
		F.init( one, 1.0 );
		typename Field::Element * Lcurr,* Rcurr,* Bcurr;
		size_t ib,  Ldim;
		int k;
		if ( Side == FflasLeft ){
			size_t j = 0;
			while ( j<R ) {
				k = ib = Q[j];
				while ((j<R) && ( (int) Q[j] == k)  ) {k++;j++;}
				Ldim = k-ib;
				Lcurr = L + j-Ldim + ib*ldl;
				Bcurr = B + ib*ldb;
				Rcurr = Lcurr + Ldim*ldl;

				ftrsm( F, Side, FflasLower, FflasNoTrans, FflasUnit, Ldim, N, one,
				       Lcurr, ldl , Bcurr, ldb );

				fgemm( F, FflasNoTrans, FflasNoTrans, M-k, N, Ldim, Mone,
				       Rcurr , ldl, Bcurr, ldb, one, Bcurr+Ldim*ldb, ldb);
			}
		}
		else{ // Side == FflasRight
			int j=R-1;
			while ( j >= 0 ) {
				k = ib = Q[j];
				while ( (j >= 0) &&  ( (int)Q[j] == k)  ) {--k;--j;}
				Ldim = ib-k;
				Lcurr = L + j+1 + (k+1)*ldl;
				Bcurr = B + ib+1;
				Rcurr = Lcurr + Ldim*ldl;

				fgemm (F, FflasNoTrans, FflasNoTrans, M,  Ldim, N-ib-1, Mone,
				       Bcurr, ldb, Rcurr, ldl,  one, Bcurr-Ldim, ldb);

				ftrsm (F, Side, FflasLower, FflasNoTrans, FflasUnit, M, Ldim, one,
				       Lcurr, ldl , Bcurr-Ldim, ldb );
			}
		}
	}


	template<class Field>
	static void trinv_left( const Field& F, const size_t N, const typename Field::Element * L, const size_t ldl,
				typename Field::Element * X, const size_t ldx ){
		for (size_t i=0; i<N; ++i)
			fcopy (F, N, X+i*ldx, 1, L+i*ldl, 1);
		ftrtri (F, FflasLower, FflasUnit, N, X, ldx);
		//invL(F,N,L,ldl,X,ldx);
	}

	template <class Field>
	static size_t KrylovElim( const Field& F, const size_t M, const size_t N,
				  typename Field::Element * A, const size_t lda, size_t*P,
				  size_t *Q, const size_t deg, size_t *iterates, size_t * inviterates, const size_t maxit,size_t virt);

	template <class Field>
	static size_t  SpecRankProfile (const Field& F, const size_t M, const size_t N,
					typename Field::Element * A, const size_t lda, const size_t deg, size_t *rankProfile);
	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	CharpolyArithProg (const Field& F, std::list<Polynomial>& frobeniusForm,
			   const size_t N, typename Field::Element * A, const size_t lda, const size_t c);

	template <class Field>
	static void CompressRows (Field& F, const size_t M,
				  typename Field::Element * A, const size_t lda,
				  typename Field::Element * tmp, const size_t ldtmp,
				  const size_t * d, const size_t nb_blocs);

	template <class Field>
	static void CompressRowsQK (Field& F, const size_t M,
				  typename Field::Element * A, const size_t lda,
				  typename Field::Element * tmp, const size_t ldtmp,
				  const size_t * d,const size_t deg, const size_t nb_blocs);

	template <class Field>
	static void DeCompressRows (Field& F, const size_t M, const size_t N,
					    typename Field::Element * A, const size_t lda,
					    typename Field::Element * tmp, const size_t ldtmp,
					    const size_t * d, const size_t nb_blocs);
	template <class Field>
	static void DeCompressRowsQK (Field& F, const size_t M, const size_t N,
					    typename Field::Element * A, const size_t lda,
					    typename Field::Element * tmp, const size_t ldtmp,
					    const size_t * d, const size_t deg, const size_t nb_blocs);

	template <class Field>
	static void CompressRowsQA (Field& F, const size_t M,
					    typename Field::Element * A, const size_t lda,
					    typename Field::Element * tmp, const size_t ldtmp,
					    const size_t * d, const size_t nb_blocs);
	template <class Field>
	static void DeCompressRowsQA (Field& F, const size_t M, const size_t N,
					      typename Field::Element * A, const size_t lda,
					      typename Field::Element * tmp, const size_t ldtmp,
					      const size_t * d, const size_t nb_blocs);


protected:


	// Subroutine for Keller-Gehrig charpoly algorithm
	// Compute the new d after a LSP ( d[i] can be zero )
	template<class Field>
	static size_t
	newD( const Field& F, size_t * d, bool& KeepOn,
	      const size_t l, const size_t N,
	      typename Field::Element * X,
	      const size_t* Q,
	      std::vector<std::vector<typename Field::Element> >& minpt);

	template<class Field>
	static size_t
	updateD(const Field& F, size_t * d, size_t k,
		std::vector<std::vector<typename Field::Element> >& minpt );

	//---------------------------------------------------------------------
	// RectangleCopyTURBO: Copy A to T, with respect to the row permutation
	//                     defined by the lsp factorization of located in
	//                     A-dist2pivot
	//---------------------------------------------------------------------
	template <class Field>
	static void
	RectangleCopyTURBO( const Field& F, const size_t M, const size_t N,
		       const size_t dist2pivot, const size_t rank,
		       typename Field::Element * T, const size_t ldt,
		       const typename Field::Element * A, const size_t lda ){

		const typename Field::Element * Ai = A;
		typename Field::Element * T1i = T, T2i = T + rank*ldt;
		size_t x = dist2pivot;
		for (; Ai<A+M*lda; Ai+=lda){
			while ( F.isZero(*(Ai-x)) ) { // test if the pivot is 0
				fcopy( F, N, T2i, 1, Ai, 1);
				Ai += lda;
				T2i += ldt;
			}
			fcopy( F, N, T1i, 1, Ai, 1);
			T1i += ldt;
			x--;
		}
	}



	//---------------------------------------------------------------------
	// LUdivine_construct: (Specialisation of LUdivine)
	// LUP factorisation of X, the Krylov base matrix of A^t and v, in A.
	// X contains the nRowX first vectors v, vA, .., vA^{nRowX-1}
	// A contains the LUP factorisation of the nUsedRowX first row of X.
	// When all rows of X have been factorized in A, and rank is full,
	// then X is updated by the following scheme: X <= ( X; X.B ), where
	// B = A^2^i.
	// This enables to make use of Matrix multiplication, and stop computing
	// Krylov vector, when the rank is not longer full.
	// P is the permutation matrix stored in an array of indexes
	//---------------------------------------------------------------------

	template <class Field>
	static size_t
	LUdivine_construct( const Field& F, const FFLAS_DIAG Diag,
			    const size_t M, const size_t N,
			    const typename Field::Element * A, const size_t lda,
			    typename Field::Element * X, const size_t ldx,
			    typename Field::Element * u, size_t* P,
			    bool computeX, const FFPACK_MINPOLY_TAG MinTag,
			    const size_t kg_mc, const size_t kg_mb, const size_t kg_j );

	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	KellerGehrig( const Field& F, std::list<Polynomial>& charp, const size_t N,
		      const typename Field::Element * A, const size_t lda );

	template <class Field, class Polynomial>
	static int
	KGFast ( const Field& F, std::list<Polynomial>& charp, const size_t N,
		 typename Field::Element * A, const size_t lda,
		 size_t * kg_mc, size_t* kg_mb, size_t* kg_j );

	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	KGFast_generalized (const Field& F, std::list<Polynomial>& charp,
			    const size_t N,
			    typename Field::Element * A, const size_t lda);


	template<class Field>
	static void
	fgemv_kgf( const Field& F,  const size_t N,
		   const typename Field::Element * A, const size_t lda,
		   const typename Field::Element * X, const size_t incX,
		   typename Field::Element * Y, const size_t incY,
		   const size_t kg_mc, const size_t kg_mb, const size_t kg_j );

	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	LUKrylov( const Field& F, std::list<Polynomial>& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  typename Field::Element * U, const size_t ldu);

	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	Danilevski (const Field& F, std::list<Polynomial>& charp,
		    const size_t N, typename Field::Element * A, const size_t lda);

	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	LUKrylov_KGFast( const Field& F, std::list<Polynomial>& charp, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * X, const size_t ldx);
};

#include "ffpack_ludivine.inl"
#include "ffpack_minpoly.inl"
#include "ffpack_charpoly_kglu.inl"
#include "ffpack_charpoly_kgfast.inl"
#include "ffpack_charpoly_kgfastgeneralized.inl"
#include "ffpack_charpoly_danilevski.inl"
#include "ffpack_charpoly.inl"
#include "ffpack_krylovelim.inl"
#include "ffpack_frobenius.inl"
#ifdef _LINBOX_LINBOX_CONFIG_H
}
#endif
#endif // __FFPACK_H
