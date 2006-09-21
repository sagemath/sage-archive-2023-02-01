freeze;

// do not freeze;
/*****************************************************************

 shark.m - Computation of upper bounds for
 		 the p-primary part of the
		 Tate-
	SHA  farevich group and
 		 the
	RanK of an elliptic curve.


 cw jul05

 ****************************************************************/

/*
	to be loaded with AttachSpec("shark.spec");
*/

/*****************************************************************
	How do you use this package?

	rankupperbound(E) returns an upper bound for the
		rank of the Mordell-Weil group of an elliptic curve
		E defined over Q using p-adic L-functions for different
		small nice primes p.

	shaupperbound(E,p) returns an upper bound for the
		size of the p-primary part of the Tate-Shafarevich group
		of the elliptic curve E defined over Q.
		In particular it proves its finiteness.
		E has to have good ordinary reduction at p > 2.

	SetVerbose("Shark",1) gives the details of the computations.

*****************************************************************/

/*****************************************************************
	Examples :

	(see also examples.m)

 * a curve with non-trivial Sha[2].
   rankupperbound can decide that the rank has to be 1.

> E:=EllipticCurve([1,-1,1,-8587,-304111]);
> MordellWeilRankBounds(E);
1 3
> rankupperbound(E);
1


 * a curve of rank 3.
   shaupperbound can prove that the 5-primary part of Sha is trivial.

> E:=EllipticCurve([ 0, 0, 1, -7, 6 ]);
> MordellWeilRank(E);
3
 > shaupperbound(E,5);
1


*****************************************************************/

/*****************************************************************
 	Problems :

	- does not work for  p = 3 when rank > 0.
	- too slow.

	to do :
	- Isogenous curves.
	- upper bound when rho is non-split Cartan or S_4.
	- p = 3 and p = 2.
	- supersingular and bad primes.

*****************************************************************/

declare verbose Shark, 1;

import "myl.m": PnmodT,cyclotomicmodT, formalpha;

// returns [r,v] where r is the first coefficient of the L-function
// that is known to be non-zero with the given precp and
// v is the valuation of the r-th coefficient.

// replaces bsdpdata in myl.m

function pinformation(V,alpha,p,precp,precT )
	R<T>:=PolynomialRing(Rationals());

	vprint Shark : "    computing the L-function.";
	P := PnmodT(V,p,precp,precT,1);
	Q := PnmodT(V,p,precp-1,precT,1);
	L := 1/alpha^precp*P-1/alpha^(precp+1)*Q*cyclotomicmodT(p,precp,precT) mod T^precT;
	vprint Shark : "    ...done.";

	L := L + 0 * T;
	if L eq 0 then
		vprint Shark : "    precisision is not enough !!";
		return([-1,-1]);
 	else
		c := Coefficient(L,0);
		if c ne 0 then
			v := Valuation(c,p);
		else
			v := 10000;
		end if;
		r := 0;
		 // check if the r-th coefficent can be distinguished from 0
		 // note that D^r L_p (0) is congruent to D^r L(0) mod p^(precp-1)
		while ((c eq 0 ) or (v ge precp -1 - Valuation(Factorial(r),p) ) ) and ( r le Degree(L) ) do
			r := r+1;
			c := Coefficient(L,r);
			v := Valuation(c,p);
		end while;
		if r gt Degree(L) then
			vprint Shark : "    precisision is not enough.";
			return([-1,-1]);
 	    else
			return([r,v]);
		end if;
	end if;

end function;


intrinsic rankupperbound(E : precT := 10) -> RngIntElt
{
 gives an upper bound for the rank of the Mordell-Weil group of E(Q).
 this is useful when the 2-descent can not conclude due to the presence of
 2-torsion in Sha.
}
/* to do :

 - using first AnalyticRank to determine when to stop ?
 - more primes ?
 - p = 2 ?

*/

	N := Conductor(E);
	vprint Shark : "Computing the spaces of modular symbols.";
	V := MSpaces(E);

	p := 3;
	rb:= 10000;
	while p lt 30 do
		if (N mod p ne 0) then
			apE := ap(E,p);
			if (apE mod p ne 0 ) then
				if is_rho_surjective(E,p) then
					// we may use the good ordinary prime p
					// feasable precp :: ??
					case p :
						when 3 : 	precp := 6;
						when 5 : 	precp := 4;
						when 7,11: 	precp := 3;
						else  		precp := 2;
					end case;
					vprint Shark : "Computing upper bound using p = ",p;
					alpha := formalpha(N,p,apE,1,11);
					pinf := pinformation(V,alpha,p,precp,precT);
					vprint Shark : "Upper bound : ", pinf[1];
					vprint Shark : "Valuation of the leading coefficient :", pinf[2];
					if (pinf[1] ne -1) and (pinf[1] lt rb) then
						rb := pinf[1];
					end if;
					// we can not get lower than that.
					if rb eq 0 then
						break;
					end if;
				end if;
			end if;
		end if;
		p:=NextPrime(p);
	end while;
	if rb eq 10000 then
		print "Could not determine an upper bound !!";
	end if;
	return(rb);
end intrinsic;




intrinsic shaupperbound(E,p) -> RngIntElt
{
 proves the finiteness of the p-primary part
 of the Tate-Shafarevich group and
 gives an upper bound for its size.
 Currently only for good ordinary primes > 2 with surjective mod p rep.
}

	require IsPrime(p) : "Argument 2 must be a prime.";
	require p gt 2 : "Argument 2 has to be an odd prime.";
	N := Conductor(E);
	require N mod p ne 0 :"Arguement 2 has to be a prime of good reduction.";
	apE := ap(E,p);
	require apE mod p ne 0 :"Argument 2 has to be an ordinary prime.";
	require is_rho_surjective(E,p) : "The mod p representation might not be surjective.";


	vprint Shark : "Computing the compex analytic information";
	ran, lstar := AnalyticRank(E);
	vprint Shark : "  analytic rank : ",ran;

// trying to compute a basis of the Mordell-Weil group.
 	if ran gt 0 then
		vprint Shark : "Computing bounds for the rank using a 2-descent.";
		a,b := MordellWeilRankBounds(E);

		M,f := MordellWeilGroup(E);
		B := [f(M.i) : i in [1..Ngens(M)] | Order(M.i) eq 0];
		vprint Shark : "  number of generators : ",a;
		vprint Shark : "  upper bound for the rank from 2-descent :",b;

		if a ne ran then
			print "Probably we did not find the full Mordell-Weil group.";
			// shall we abort ?
		end if;
	else
		a:=0;
		b:=0;
	end if;

// compute the expected order of sha(p)
	if a gt 0 then
		reginf := Regulator(B);
	else
		reginf := 1;
	end if;
	shan := lstar * (#TorsionSubgroup(E))^2/reginf/RealPeriod(E)/&*[TamagawaNumber(E,q) : q in BadPrimes(E)];
	if Discriminant(E) gt 0 then
		shan := shan/2;
	end if;
	shan := Round(shan);
	shan := p^Valuation(shan,p);
	vprint Shark : "  the analytic order of the ",p,"-primary part of Sha : ",shan;

//  later we can short-cut here for non-vanishing L-functions,
// once the algorithm is stable

// compute information from the p-adic L- function.

	vprint Shark : "Computing the analytic p-adic information";
	alpha := formalpha(N,p,apE,1,11);
	vprint Shark : "  computing the spaces of modular symbols.";
	V := MSpaces(E);
	pinf := [-1,-1];
	// this estimate should be sharp if the p-adic regulator
	// is a unit in Z_p (or p^-2 Z_p* if anomalous)
	// otherwise we have to increase the
	// precp in order to distinguish the leading term from 0
	precp := 1 + Valuation(BestApproximation(Factorial(ran) * lstar /RealPeriod(E)/reginf,100),p);
	while (pinf[1] eq -1 ) and precp lt 7 do
	    precp := precp + 1;
		vprint Shark : "  precp increased to ",precp;
		pinf := pinformation(V,alpha,p,precp,ran+1);
	end while;
	vprint Shark : "  order of vanishing of the ",p,"-adic L-function : ",pinf[1];
	vprint Shark : "  valuation of the leading term : ",pinf[2];

	if a eq b then
		vprint Shark : "Rank is equal to", a;
		r := a;
	elif pinf[1] eq a then
		vprint Shark : "Rank is equal to", a;
		r := a;
	else
		b2 := rankupperbound(E : precT := ran+1);
		if b2 eq a then
			vprint Shark : "Rank is equal to", a;
			r := a;
		else
			vprint Shark : "Could not determine the rank. Abort.";
			r := -1;
				// for rank 1 we could add the use of Heegner points ??
		end if;
	end if;

    if r eq -1 then
		return(-1); // returns -1 if we can not compute a basis of E(Q).
	else

			// do we have enough prec in pinf to compute ?
		vprint Shark : "Is the precp high enough ?";
		oldprecp := precp;
		while pinf[1] ne r do
		   	precp := precp + 1;
			vprint Shark : "  precp increased to ",precp;
			pinf := pinformation(V,alpha,p,precp,r+1);
		end while;
		if precp gt oldprecp then
			vprint Shark : "  order of vanishing of the ",p,"-adic L-function : ",pinf[1];
			vprint Shark : "  valuation of the leading term : ",pinf[2];
			vprint Shark : "  the ",p,"-priary part of Sha is finite !!";
		else
			vprint Shark : "  yes.";
		end if;

		// now we should have pinf = [r, val of leading term = rth term]


// computation of the p-adic regulator
		vprint Shark : "Computing the p-adic regulator.";

		if r gt 0 then
			if p eq 3 then
				h := height_function(E, p, Max(precp+2,10) : use_integrality := true );
				// does not seem to work either
			else
				h := height_function(E, p, Max(precp+2,10) );
			end if;
			K := Parent(h(B[1]));
			M := MatrixAlgebra(K,#B)!0;
			function pairing(P,Q)
				return (h(P+Q)-h(P)-h(Q))/2;
			end function;
			for i in [1..#B] do
				for j in [i..#B] do
					aij := pairing(B[i],B[j]);
					M[i,j] := aij;
					M[j,i] := aij;
				end for;
			end for;
			// we assume silently that the reg is non-zero.
			// this must follow from Kato-... and the above computation
			// that pinf[1] = r
			regp := Valuation(Determinant(M));
	 		vprint Shark : "  the normalized regulator has valuation : ",regp;
		else
			regp := 0;
	 		vprint Shark : "  the normalized regulator has valuation : ",regp;
		end if;

// putting things together
		 Np := apE -p -1;

		vprint Shark : "Putting things together : ";
		vprint Shark : "  by Kato- Perrin-Riou/Schneider we have";
		vprint Shark : "  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (# E(Q)_tors)";
		vprint Shark : "  ",pinf[2]," >= vp(# Sha) + ",Valuation(&*[TamagawaNumber(E,q) : q in BadPrimes(E)],p)," + 2 *",Valuation(Np,p)," + ",regp," - 2 * ",Valuation(#TorsionSubgroup(E),p);


		 shap := pinf[2] + 2* Valuation(#TorsionSubgroup(E),p);
		 shap := shap - 2* Valuation(Np,p)- Valuation(&*[TamagawaNumber(E,q) : q in BadPrimes(E)],p);
		 shap := shap - regp;
		 shap := p^shap;

		vprint Shark : "  so #Sha(p) <= ", shap;
		if shap eq 1 then
			vprint Shark : "  hence the main conjecture holds.";
		else
			vprint Shark : "  this is an equality under the assumption of the main conjecture (...Urban-Skinner).";
		end if;

		// check.
		if shap ne shan then
			vprint Shark : "  but there is a contradiction with B-SwD !!";
		end if;

		return(shap);

	 end if;

end intrinsic;


