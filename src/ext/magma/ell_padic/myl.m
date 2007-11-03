freeze;

/*****************************************************************

 myl.m - Computation of the p-adic L-series
 		 for good ordinary primes only !!
		 and p > 2 only !!

 (part of shark)

 copied and adjusted from Robert Pollack's code

 cw july05

 ****************************************************************/

/**** to do

- we should have a p-adic evaluation of the modular symbol so that we could
  force the prec of it to be low

- everything for p = 2

- for multiplicative reduction

- for supersingular primes

- for additive reduction ??

****/

R<T>:=PolynomialRing(Rationals());

intrinsic ap(E::CrvEll, p::RngIntElt) ->RngIntElt
{
  Returns #E(F_p) - (p + 1).
}
   return TraceOfFrobenius(ChangeRing(E,GF(p)));
end intrinsic;

intrinsic ModularIntegral(M,r,D) -> FldRatElt
{
 Integrates the modular form attached to Argument 1 twisted by Argument 3 from Argument 2 to I * infinity
}


	function MSEven(M,r)
 		return(ModularSymbolEven(M,r)[1]/4);
	end function;

	function MSOdd(M,r)
	 	return(ModularSymbolOdd(M,r)[1]/4);
	end function;

	answer:=0;
	for a in [0..Abs(D)-1] do
		if KroneckerSymbol(D,a) ne 0 then
			if D ge 0 then
				ms:=MSEven(M[1],r-a/D);
			else
				ms:=MSOdd(M[2],r+a/D);
			end if;
			answer:=answer+ms*KroneckerSymbol(D,a);
		end if;
	end for;

	return(answer);
end intrinsic;

intrinsic MSpaces(E) ->.
{
  Return the spaces of plus and minus Modularsymbols attached to E
}

	N    := Conductor(E);
	R<x> := PolynomialRing(Rationals());
	M1   := ModularSymbols(N,2,+1);
	p:=2;
	while N mod p eq 0 do
		p:=NextPrime(p);
	end while;
	while Dimension(M1) gt 1 do
		M1:=Kernel([<p,x-ap(E,p)>],M1);
		p:=NextPrime(p);
		while N mod p eq 0 do
			p:=NextPrime(p);
		end while;
	end while;
	M2   := ModularSymbols(N,2,-1);
	p:=2;
	while N mod p eq 0 do
		p:=NextPrime(p);
	end while;
	while Dimension(M2) gt 1 do
		M2:=Kernel([<p,x-ap(E,p)>],M2);
		p:=NextPrime(p);
		while N mod p eq 0 do
			p:=NextPrime(p);
		end while;
	end while;

	return([M1,M2]);

end intrinsic;

// computes the root of x^2 -ap x +p in Q_p
function formalpha(N,p,a_p,D,m) //->RngIntElt

	vprint Shark : "    computing the root of Frobenius...";

	assert N mod p^2 ne 0;
	assert a_p mod p ne 0;
	assert D mod p ne 0;

	if N mod p eq 0 then
		temp := a_p*KroneckerSymbol(D,p);
	else
		a_p := a_p * KroneckerSymbol(D,p);
		temp := a_p mod p;
		for k in [2..m] do
			t := (-(temp^2-a_p*temp+p) div p^(k-1) mod p) * InverseMod(2*temp - a_p,p);
			temp := temp + t*p^(k-1);
		end for;
	end if;
	vprint Shark : "    ...done. alpha =",temp;


	return(temp);
end function;

intrinsic teichmuller(a,p,m)->RngIntElt
{
 Returns the Teichmuller representative of a mod p^m
}
	require IsPrime(p) : "Argument 2 must be prime.";
	if p ne 2 then
		r := a mod p;     // newton algorithm for approximating X^(p-1) - 1 == 0
		for k in [1..m] do
			r := r - (Modexp(r,p-1,p^m) - 1) * Modinv( p-1,  p^m )* Modexp(r,2-p,p^m);
			r := r mod p^m;
		end for;
	else
		if a mod 4 eq 1 then
			return(1);
		else
			return(-1);
		end if;
	end if;

	return(r);

end intrinsic;

function listofteichmuller(p,m)
	local answer;

	answer := [];
		for a in [1..p-1] do
			Append(~answer,teichmuller(a,p,m));
		end for;
	return(answer);
end function;


function cyclotomic(p,n)
	if n gt 0 then
		answer:=0;
		for a in [0..p-1] do
			answer:=answer+(1+T)^(a*p^(n-1));
		end for;
		return(answer);
	else
		return(T);
	end if;
end function;

function Pn(V,p,n,D)
	lteich := listofteichmuller(p,n+1);
	answer:=0;
	oneplusTpower:=1;
	gammapower:=1;
	gamma:=1+p;
	if p ne 3 then
		for j in [0..p^n-1] do
			for a in [1..((p-1) div 2)] do
				ms := ModularIntegral(V,gammapower * lteich[a]/p^(n+1),D);
				answer := (answer + ms * oneplusTpower) ;
			end for;
			gammapower := gammapower * gamma;
			oneplusTpower := (oneplusTpower * (1+T));
		end for;
	else
		for j in [0..p^n-1] do
			ms := ModularIntegral(V,gammapower/p^(n+1),D);
			answer := answer + ms * oneplusTpower;
			gammapower := gammapower * gamma;
			oneplusTpower := oneplusTpower * (1+T);
		end for;
	end if;

	return(answer);
end function;


intrinsic Lpfunction(E,p,n : twistD:=1) -> .
{
 Returns an approximation to the p-adic L-series.
}
	require IsPrime(p) : "Argument 2 must be prime" ;
	require p gt 2 : "You might burn in hell for doing Iwasawa theory for p = 2" ;
	require n ge 1 : "Argument 3 must be positive." ;

	R<T>:=PolynomialRing(Rationals());
	N := Conductor(E);
	assert N mod p ne 0;
	apE := ap(E,p);
	assert apE mod p ne 0;
	alpha := formalpha(N,p,apE,twistD,11);
	V := MSpaces(E);
	P := Pn(V,p,n,twistD);
	Q := Pn(V,p,n-1,twistD);
	L := 1/alpha^n*P-1/alpha^(n+1)*Q*cyclotomic(p,n);

	return L;
end intrinsic;

/**************** the same but all computations mod T^some power ****************************/


function cyclotomicmodT(p,n,precT)
	if n gt 0 then
		answer:=0;
		for a in [0..p-1] do
			answer:=answer+(1+T)^(a*p^(n-1)) mod T^precT;
		end for;
		return(answer);
	else
		return(T);
	end if;
end function;

function PnmodT(V,p,n,precT,D)
	lteich := listofteichmuller(p,n+1);
	answer:=0;
	oneplusTpower:=1;
	gammapower:=1;
	gamma:=1+p;
	if p ne 3 then
		for j in [0..p^n-1] do
			for a in [1..((p-1) div 2)] do
				ms := ModularIntegral(V,gammapower * lteich[a]/p^(n+1),D);
				answer := (answer + ms * oneplusTpower) mod T^precT ;
			end for;
			gammapower := gammapower * gamma;
			oneplusTpower := (oneplusTpower * (1+T)) mod T^precT;
		end for;
	else
		for j in [0..p^n-1] do
			ms := ModularIntegral(V,gammapower/p^(n+1),D);
			answer := answer + ms * oneplusTpower mod T^precT ;
			gammapower := gammapower * gamma;
			oneplusTpower := oneplusTpower * (1+T) mod T^precT ;
		end for;
	end if;

	return(answer);
end function;

intrinsic LpfunctionmodT(E,p,n,precT : twistD:=1) -> .
{
 Returns an approximation to the p-adic L-series.
}
	require IsPrime(p) : "Argument 2 must be prime" ;
	require p gt 2 : "You might burn in hell for doing Iwasawa theory for p = 2" ;
	require n ge 1 : "Argument 3 must be positive." ;

	R<T>:=PolynomialRing(Rationals());
	N := Conductor(E);
	assert N mod p ne 0;
	apE := ap(E,p);
	assert apE mod p ne 0;
	alpha := formalpha(N,p,apE,twistD,11);
	V := MSpaces(E);
	vprint Shark : "    computing the L-function.";
	P := PnmodT(V,p,n,precT,twistD);
	Q := PnmodT(V,p,n-1,precT,twistD);
	L := 1/alpha^n*P-1/alpha^(n+1)*Q*cyclotomicmodT(p,n,precT) mod T^precT;

	return L;
end intrinsic;



intrinsic bsdpdata(E,p,precp,precT) -> .
{
 Returns an upper bound b for the rank of E
 and the p-adic valuation of the b-th coefficient of L_p(T).
 }
// Computations are done only up to T^precT and p^precp

	require IsPrime(p):"Argument 2 must be prime" ;
	require p gt 2 : "You might burn in hell for doing Iwasawa theory for p = 2" ;
	require Conductor(E) mod p ne 0 :"The elliptic curve must have good reduction at p.";
	apE := ap(E,p);
	require apE mod p ne 0 :"The elliptic curve must have ordinary reduction at p.";


    L :=  LpfunctionmodT(E,p,precp,precT);
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

end intrinsic;


