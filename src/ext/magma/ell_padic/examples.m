freeze;

 /* EXAMXPLES FOR SHARK.M ***************************

 this is a list of interesting examples for the computations
 done by shark.m

 (the detailed output does not always correspond to the latest
 version of the print commands...)


*****************************************************/


 DB := CremonaDatabase();
 SetVerbose("Shark",1);


 /*************** some standard curves ***************/


 FirstCurve := EllipticCurve(DB,"11A1");

/* > time rankupperbound(FirstCurve);
Computing the spaces of modular symbols.
Computing upper bound using p =  3
    computing the root of Frobenius...
    ...done. alpha = 47612
    computing the L-function.
    ...done.
Upper bound :  0
Valuation of the leading coefficient : 0
0
Time: 2.850
*/

/* > time shaupperbound(FirstCurve,3);
Computing the compex analytic information
  analytic rank :  0
  the analytic order of the  3 -primary part of Sha :  1
Computing the analytic p-adic information
    computing the root of Frobenius...
    ...done. alpha = 47612
  computing the spaces of modular symbols.
  precp increased to  2
    computing the L-function.
    ...done.
  order of vanishing of the  3 -adic L-function :  0
  valuation of the leading term :  0
Rank is equal to 0
Is the precp high enough ?
  yes.
Computing the p-adic regulator.
  the normalized regulator has valuation :  0
Putting things together :
  by Kato- Perrin-Riou/Schneider we have
  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (#
E(Q)_tors)
   0  >= vp(# Sha) +  0  + 2 * 0  +  0  - 2 *  0
  so #Sha(p) <=  1
  hence the main conjecture holds.
1
Time: 0.060
*/

/*> time shaupperbound(FirstCurve,7);
1
Time: 0.460
> time shaupperbound(FirstCurve,13);
1
Time: 3.330
> time shaupperbound(FirstCurve,17);
1
Time: 7.810
> time shaupperbound(FirstCurve,23);
1
Time: 20.670
*/

 RankOneCurve := EllipticCurve(DB,"37A");

/* 13 is irregular,
	so the 2nd approximation to the rank is not sufficient.
*/

/*> time shaupperbound(FirstCurve,23);
1
Time: 20.670
> SetVerbose("Shark",1);
> rankupperbound(RankOneCurve);
Computing the spaces of modular symbols.
Computing upper bound using p =  5
    computing the root of Frobenius...
    ...done. alpha = 51341613
    computing the L-function.
    ...done.
Upper bound :  1
Valuation of the leading coefficient : 0
Computing upper bound using p =  7
    computing the root of Frobenius...
    ...done. alpha = 125426580
    computing the L-function.
    ...done.
Upper bound :  1
Valuation of the leading coefficient : 0
Computing upper bound using p =  11
    computing the root of Frobenius...
    ...done. alpha = 78049571232
    computing the L-function.
    ...done.
Upper bound :  1
Valuation of the leading coefficient : 0
Computing upper bound using p =  13
    computing the root of Frobenius...
    ...done. alpha = 2666354398667
    computing the L-function.
    ...done.
Upper bound :  3
Valuation of the leading coefficient : 0
Computing upper bound using p =  23
    computing the root of Frobenius...
    ...done. alpha = 297453551161082
    computing the L-function.
    ...done.
Upper bound :  1
Valuation of the leading coefficient : 0
Computing upper bound using p =  29
    computing the root of Frobenius...
    ...done. alpha = 56214124488451581
    computing the L-function.
    ...done.
Upper bound :  1
Valuation of the leading coefficient : 0
1
*/

/*> shaupperbound(RankOneCurve,5);
1
> shaupperbound(RankOneCurve,7);
1
> shaupperbound(RankOneCurve,11);
1
> time shaupperbound(RankOneCurve,13);
1
Time: 55.390
*/

 RankOneCurve2 := EllipticCurve(DB, "43A");

/*> shaupperbound(RankOneCurve2,5);
1
> shaupperbound(RankOneCurve2,11);
1
> shaupperbound(RankOneCurve2,13);
1
> shaupperbound(RankOneCurve2,17);
1
> shaupperbound(RankOneCurve2,19);
1
> shaupperbound(RankOneCurve2,23);
1
*/

 RankTwoCurve := EllipticCurve(DB,"389A");

/* > shaupperbound(E,5);
Computing the compex analytic information
  analytic rank :  2
Computing bounds for the rank using a 2-descent.
  number of generators :  2
  upper bound for the rank from 2-descent : 2
  the analytic order of the  5 -primary part of Sha :  1
Computing the analytic p-adic information
    computing the root of Frobenius...
    ...done. alpha = 2593532
  computing the spaces of modular symbols.
  precp increased to  2
    computing the L-function.
    ...done.
  order of vanishing of the  5 -adic L-function :  2
  valuation of the leading term :  0
Rank is equal to 2
Is the precp high enough ?
  yes.
Computing the p-adic regulator.
  the normalized regulator has valuation :  0
Putting things together :
  by Kato- Perrin-Riou/Schneider we have
  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (#E(Q)_tors)
   0  >= vp(# Sha) +  0  + 2 * 0  +  0  - 2 *  0
  so #Sha(p) <=  1
  hence the main conjecture holds.
1
*/

 RankThreeCurve := EllipticCurve(DB, "5077A");

/*> shaupperbound(E,5);
Computing the compex analytic information
  analytic rank :  3
Computing bounds for the rank using a 2-descent.
  number of generators :  3
  upper bound for the rank from 2-descent : 3
  the analytic order of the  5 -primary part of Sha :  1
Computing the analytic p-adic information
    computing the root of Frobenius...
    ...done. alpha = 23157316
  computing the spaces of modular symbols.
  precp increased to  2
    computing the L-function.
    ...done.
  order of vanishing of the  5 -adic L-function :  3
  valuation of the leading term :  0
Rank is equal to 3
Is the precp high enough ?
  yes.
Computing the p-adic regulator.
  the normalized regulator has valuation :  -2
Putting things together :
  by Kato- Perrin-Riou/Schneider we have
  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (#E(Q)_tors)
   0  >= vp(# Sha) +  0  + 2 * 1  +  -2  - 2 *  0
  so #Sha(p) <=  1
  hence the main conjecture holds.
1
*/

 MyFavouriteCurve := EllipticCurve(DB, "5692A");
	// rank 2 with L_3 = 3^3 * T^2 + ...
	// should give trivial Sha(3).
	// but 3-adic heights are not implemented :
	// solve_for_s2(
	//    sigma: t + (s2 + 1/3)*t^3 + (1/2*s2^2 + s2 - 263/36)*t^5 + (1/6*s2^...,
	//    p: 3
	//)
	//In file "/home/chrigu/shark/padic_height.m", line 383, column 10:
	//>>          assert #X le 1;
	//            ^
	//Runtime error in assert: Assertion failed


/***************  Examples from my fine Sha preprint... ********/

 FSha1 := EllipticCurve(DB, "930O3");
 FSha2 := EllipticCurve(DB, "210E5");
 FSha3 := EllipticCurve(DB, "1287E2");
	 //> MordellWeilRankBounds(FSha1);
	 //0 2

	//> time rankupperbound(FSha1);
	//Computing the spaces of modular symbols.
	//Computing upper bound using p =  11
	//    computing the root of Frobenius...
	//    ...done. alpha = 1435809126759
	//    computing the L-function.
	//    ...done.
	//Upper bound :  0
	//Valuation of the leading coefficient : 0
	//0
	//Time: 39.740

	//> rankupperbound(FSha2);
	//    computing the spaces of modular symbols.
	//    ...done.
	//Computing upper bound using p =  11
	//    computing the root of Frobenius...
	//    ...done. alpha = 1435809126759
	//    computing the L-function.
	//    ...done.
	//Upper bound :  0
	//Valuation of the leading coefficient : 0
	//0

	// > time rankupperbound(FSha3);
	//Computing the spaces of modular symbols.
	//Computing upper bound using p =  5
	//    computing the root of Frobenius...
	//    ...done. alpha = 143970887
	//    computing the L-function.
	//    ...done.
	//Upper bound :  0
	//Valuation of the leading coefficient : 0
	//0
	//Time: 10.350

 FSha4twist := EllipticCurve(DB, "17A2");
 FSha4 := MinimalModel(QuadraticTwist (FSha4twist, -3*13*157));
	 // > MordellWeilRankBounds(FSha4);
	 // 0 4
	 // due to Sha[2] = Z/2 ^4
	 // conductor too big for modularSymbols
	//	 Computing the spaces of modular symbols.

	//Current total memory usage: 143.3MB, failed memory request: 11670.2MB
	//System error: Out of memory.
	//All virtual memory has been exhausted so Magma cannot perform this statement.

 	// we might be able to this usind twistD parameter of Pollack's L-functions.


 /*****************  from allbigSha of Cremona **********************/

 // the lowest conductor example of non-trivial p-primary parts of sha
 // such that E has good ordinary reduction at p...
 // shan has order n^2.

 Sha3 :=EllipticCurve( [1, 1, 0, -39, -112]);

/* > shaupperbound(Sha3,3);
Computing the compex analytic information
  analytic rank :  0
  the analytic order of the  3 -primary part of Sha :  9
Computing the analytic p-adic information
    computing the root of Frobenius...
    ...done. alpha = 2696
  computing the spaces of modular symbols.
  precp increased to  2
    computing the L-function.
    ...done.
    precisision is not enough.
  precp increased to  3
    computing the L-function.
    ...done.
    precisision is not enough.
  precp increased to  4
    computing the L-function.
    ...done.
  order of vanishing of the  3 -adic L-function :  0
  valuation of the leading term :  2
Rank is equal to 0
Is the precp high enough ?
  yes.
Computing the p-adic regulator.
  the normalized regulator has valuation :  0
Putting things together :
  by Kato- Perrin-Riou/Schneider we have
  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (#
E(Q)_tors)
   2  >= vp(# Sha) +  0  + 2 * 0  +  0  - 2 *  0
  so #Sha(p) <=  9
  this is an equality under the assumption of the main conjecture
(...Urban-Skinner).
9
*/

Sha3second :=EllipticCurve([1, 1, 0, -34, -135]);

/*> shaupperbound(Sha3second,3);
Computing the compex analytic information
  analytic rank :  0
  the analytic order of the  3 -primary part of Sha :  9
Computing the analytic p-adic information
    computing the root of Frobenius...
    ...done. alpha = 2696
  computing the spaces of modular symbols.
  precp increased to  2
    computing the L-function.
    ...done.
    precisision is not enough.
  precp increased to  3
    computing the L-function.
    ...done.
    precisision is not enough.
  precp increased to  4
    computing the L-function.
    ...done.
  order of vanishing of the  3 -adic L-function :  0
  valuation of the leading term :  2
Rank is equal to 0
Is the precp high enough ?
  yes.
Computing the p-adic regulator.
  the normalized regulator has valuation :  0
Putting things together :
  by Kato- Perrin-Riou/Schneider we have
  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (#
E(Q)_tors)
   2  >= vp(# Sha) +  0  + 2 * 0  +  0  - 2 *  0
  so #Sha(p) <=  9
  this is an equality under the assumption of the main conjecture
(...Urban-Skinner).
9
*/

 Sha9 :=EllipticCurve( [1, 0, 1, -22376989924, -1303266065395934]);

/*   ??   */

 Sha5  := EllipticCurve([1,-1,0,-332311,-73733731]);

/*> shaupperbound(Sha5,5);
Computing the compex analytic information
  analytic rank :  0
  the analytic order of the  5 -primary part of Sha :  25
Computing the analytic p-adic information
    computing the root of Frobenius...
    ...done. alpha = 143970887
  computing the spaces of modular symbols.
  precp increased to  4
    computing the L-function.
    ...done.
  order of vanishing of the  5 -adic L-function :  0
  valuation of the leading term :  2
Rank is equal to 0
Is the precp high enough ?
  yes.
Computing the p-adic regulator.
  the normalized regulator has valuation :  0
Putting things together :
  by Kato- Perrin-Riou/Schneider we have
  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (#E(Q)_tors)
   2  >= vp(# Sha) +  0  + 2 * 0  +  0  - 2 *  0
  so #Sha(p) <=  25
  this is an equality under the assumption of the main conjecture (...Urban-Skinner).
25 */

 Sha7  :=EllipticCurve([0, 0, 0, -4062871, -3152083138]);

/*> shaupperbound(Sha7,7);
Computing the compex analytic information
  analytic rank :  0
  the analytic order of the  7 -primary part of Sha :  49
Computing the analytic p-adic information
    computing the root of Frobenius...
    ...done. alpha = 3379162294
  computing the spaces of modular symbols.
  precp increased to  4
    computing the L-function.
    ...done.
  order of vanishing of the  7 -adic L-function :  0
  valuation of the leading term :  2
Rank is equal to 0
Is the precp high enough ?
  yes.
Computing the p-adic regulator.
  the normalized regulator has valuation :  0
Putting things together :
  by Kato- Perrin-Riou/Schneider we have
  vp(L_p*(0)) >= vp(# Sha) + vp(prod c_p) + 2 * vp(N_p) + vp(Reg_p) - 2 * vp (#
E(Q)_tors)
   2  >= vp(# Sha) +  0  + 2 * 0  +  0  - 2 *  0
  so #Sha(p) <=  49
  this is an equality under the assumption of the main conjecture
(...Urban-Skinner).
49
*/

 Sha11 := EllipticCurve([1,1,0,-16294012950,-800559991923500]);
 Sha13 := EllipticCurve([1,-1,1,-911138880,-10586098442003]);
 Sha17 := EllipticCurve([1,0,1,-59851665231,-5635904962581246]);

/*   ??   */

