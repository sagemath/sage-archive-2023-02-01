freeze;

/*****************************************************************

 gl2 -- algorithm to prove that the mod p representation is
        surjective

 (part of shark)

 cw june 05

 ****************************************************************/

/*****************************************************************

 to do : this algorithm can probably be made more efficient...


 ****************************************************************/

intrinsic is_rho_surjective(E::CrvEll, p::RngIntElt) -> .
{
 given an elliptic curve E and a prime p , this algorithm tries to prove that
 the mod p representation G_Q -> GL_2(F_p) is surjective
  true : if it is surjective
  false : likely that it is not surjective

}
	require IsPrime(p) : "p must be a prime";
	_ , d := IsogenousCurves(E);
	if Gcd(d,p) ne 1 then
		return false;   // contained in a Borel or a split Cartan
	else
	if p eq 2 then
		return not IsSquare(Discriminant(E));  // is smaller only when Delta is a square (serre 5.3)
	else
	if p eq 3 then
		return not IsCube(Discriminant(E)); // serre 5.3.
	else
		if IsSquarefree(Conductor(E)) then
			return true;   // serre shows that for semistable curves either there is a p-isogeny or it is surjective
		else
			boo:=false;
			if not(IsIntegral(jInvariant(E))) then
				f := Factorisation(Denominator(jInvariant(E)));
				for ff in f do
					ell := ff[1];
					if Gcd(-Valuation(jInvariant(E),ell), p) eq 1 then
						boo:=true;  // inerta at ell contains an element of order p -> surjective
					end if;
				end for;
			end if;
					// now using Prop 19 in serre
			ell := 2; booi := false; booii:= false; booiii:=false;
			while not boo and ell lt 1000 do
				if Gcd(ell, p*Conductor(E) ) eq 1 then
					aell := ap(E,ell);
					if  aell mod p ne 0 and LegendreSymbol(aell^2 - 4*ell,p) eq 1 then booi :=true ; end if;
					if  aell mod p ne 0 and LegendreSymbol(aell^2 - 4*ell,p) eq -1 then booii :=true ; end if;
					if ell ne p then
						_, ellinv, _ := Xgcd(ell,p);
						uell := aell^2*ellinv mod p;
						if uell^2-3*uell+1 mod p eq 0 and uell notin [0,1,2,4] then
							booiii:=true;
						end if;
					end if;
					boo:=booi and booii and booiii;
				end if;
				ell:=NextPrime(ell);
   			end while;
			return boo ;  // in practice this has a very week chance of not being able to overlook cases when it is indeed surjective
						// there should be a proveable bound.
		end if;
	end if;
	end if;
	end if;
end intrinsic;


intrinsic ap(E::CrvEll, p::RngIntElt) ->RngIntElt
{
  Returns #E(F_p) - (p + 1).
}
   return TraceOfFrobenius(ChangeRing(E,GF(p)));
end intrinsic;

intrinsic IsCube(n::FldRatElt) -> BoolElt
{
 Returns true is n is a cube
 !! this is not very clever but it works as long as n is an integer not too large.
}
 return IsIntegral(Abs(n)^(1/3));
end intrinsic;
