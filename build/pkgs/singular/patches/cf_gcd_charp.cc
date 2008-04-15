#include "config.h"
#include <stdio.h>
#ifndef NOSTREAMIO
#include <iostream>
#endif
#include <math.h>
#include <cf_factory.h>
#include <variable.h>
#include <canonicalform.h>
#include <cf_iter.h>
#include <cf_random.h>
#include <cf_algorithm.h>
#include <cf_map.h>

static CanonicalForm contentWRT(const CanonicalForm & F, const int lev);
static int degWRT(CanonicalForm F, const int lev);
static CanonicalForm lcoefWRT( const CanonicalForm & F, const int lev);
static CanonicalForm newtonInterp(CFList &alpha, CFList &u, const Variable & x);
static CanonicalForm simpleGCD(const CanonicalForm & A, const CanonicalForm & B);
CanonicalForm newGCD(CanonicalForm A, CanonicalForm B);
static CanonicalForm GFPowDown(const CanonicalForm & A, int k);
static CanonicalForm GFPowUp(const CanonicalForm & A, int k);
static CanonicalForm GFMapUp(const CanonicalForm & A, int k);
static CanonicalForm GFMapDown(const CanonicalForm & A, int k);




static CanonicalForm
contentWRT0(const CanonicalForm & F, const int lev, const int lev0)
// Computes the content of a polynomial, considering the variables of level
// larger than lev as as parameters.
{
  if(F.inBaseDomain() || (F.mvar().level() <= lev))
  {
    return(1);
  }
  else
  {
    if (F.isUnivariate())
    {
      return(F);
    }
    else
    {
      CanonicalForm pol;
      CFIterator i = CFIterator(F);
      CanonicalForm c = 0;
      for(; i.hasTerms(); i++)
      {
        CanonicalForm cc=i.coeff();
        if (cc.level() > lev0)
          pol = contentWRT0(cc, lev, lev0-1 );
        else
          pol = contentWRT(cc, lev -1);
        c = gcd(c, pol);
        if(c.isOne())
        {
          return c;
        }
      }
      return c;
    }
  }
}

static CanonicalForm
contentWRT(const CanonicalForm & F, const int lev)
// Computes the content of a polynomial, considering the variables of level
// larger than lev as as parameters.
{
  //if (lev <0)
  //  cout << "contentWRT:" <<F << " " << lev <<"\n";//printf("level=%d\n",lev);
  //int debug = 0;
  //if(debug == 1)cout << "contentWRT input - F:  " << F << " lev: " << lev << endl;
  if(F.inBaseDomain() || (F.mvar().level() <= lev))
  {
    //if(debug == 1)cout << "output 1:" << F << endl;
    return(1);
  }
  else
  {
    if (F.isUnivariate())
    {
      //if(debug == 1)cout << "output 2:" << F << endl;
      return(F);
    }
    else
    {
      CanonicalForm pol;
      CFIterator i = CFIterator(F);
      CanonicalForm c = 0;
      for(; i.hasTerms(); i++)
      {
        CanonicalForm cc=i.coeff();
        if (cc.level() > lev)
          pol = contentWRT0(cc, lev, cc.level());
        else
          pol = contentWRT(i.coeff(), lev - 1);
        //if(debug == 1)cout << "c: " << c << " - pol : " << pol << endl;
        c = gcd(c, pol);
        //if(debug == 1)cout << "c temp:" << c << endl;
        if(c.isOne())
        {
          //if(debug == 1)cout << "output 3:" << c << endl;
          return c;
        }
      }
      //if(debug == 1)cout << "output 4:" << c << endl;
      return c;
    }
  }
}

// Degree of F wrt all variables of level smaller than or equal to lev
static int degWRT(CanonicalForm F, const int lev)
{
  int deg = F.degree(lev);
  for(int i = lev; i >= 1; i--)
  {
    F = F.LC(i);
    deg = deg + F.degree(i - 1);
  }
  return deg;
}

// Leading coefficient of F wrt all variables of level smaller than or equal to lev.
static CanonicalForm lcoefWRT(const CanonicalForm & F, const int lev)
{
  CanonicalForm FL;
  FL = F;

  for(int i = lev; i >= 1; i--)
  {
    FL = FL.LC(i);
  }
  return FL;
}


static CanonicalForm
newtonInterp(const CanonicalForm alpha, const CanonicalForm u, const CanonicalForm newtonPoly, const CanonicalForm oldInterPoly, const Variable & x)
// Newton interpolation - Incremental algorithm.
// Given a list of values alpha_i and a list of polynomials u_i, 1 <= i <= n,
// computes the interpolation polynomial assuming that
// the polynomials in u are the results of evaluating the variabe x
// of the unknown polynomial at the alpha values.
// This incremental version receives only the values of alpha_n and u_n and
// the previous interpolation polynomial for points 1 <= i <= n-1, and computes
// the polynomial interpolating in all the points.
// newtonPoly must be equal to (x - alpha_1) * ... * (x - alpha_{n-1})
{
  /*
  int debug = 0;
  if(debug)
  {
    cout << "newtonInterpIncremental input - variable " << x <<endl;
    cout << "newtonInterpIncremental input - alpha:" << alpha << " - u: " << u << endl;
    cout << "newtonInterpIncremental input - newtonPoly: " << newtonPoly << " - oldInterPoly: " << oldInterPoly << endl;
  }
  */

  CanonicalForm interPoly;

  interPoly = oldInterPoly + (u - oldInterPoly(alpha, x)) / newtonPoly(alpha, x) * newtonPoly;

  //if(debug)cout << "newtonInterpIncremental output:" << interPoly << endl;
  return interPoly;
}


// Univariate GCD
// Naive implementation of Euclidean Algorithm.
// Should be used only for GCD of univariate polys over finite fields,
// where there is no coefficient growth.
static CanonicalForm simpleGCD(const CanonicalForm & A, const CanonicalForm & B)
{
  //int debug = 0;
  CanonicalForm P1, P2;
  /*
  if(debug == 1)
  {
    cout << "simpleGCD input" << endl;
    cout << "A: " << A << endl;
    cout << "B: " << B << endl;
  }
  */
  if (A.degree() > B.degree())
  {
    P1 = A;
    P2 = B;
  }
  else
  {
    P1 = B;
    P2 = A;
  }
  CanonicalForm rem;
  while(!P2.isZero())
  {
    // cout << "P1: " << P1 << " - P2: " << P2 << endl;
    rem = P1 % P2;
    // cout << "rem: " << rem << endl;
    P1 = P2;
    P2 = rem;
  }
  //if(debug == 1)cout << "simpleGCD output: " << P1 << endl;
  return(P1);
}


static CanonicalForm GFPowDown(const CanonicalForm & A, int k)
// Replaces all coefficients of A in the ground field by their k-th roots.
// It assumes that the k-th roots exist.
// This procedure is used to map down a polynomial from an extension field.
{
  CanonicalForm result = 0;
  int i, j;
  int fieldSize = pow(getCharacteristic(), getGFDegree());
  CanonicalForm g;
  for(i = 0; i <= degree(A); i++)
  {
    if(!A[i].isZero())
    {
      //cout << "i: " << i << endl;
      //cout << "A[i]: " << A[i] << endl;
      if(A[i].inBaseDomain())
      {
        //cout << "GFPowDown - inBaseDomain" << endl;
        g = getGFGenerator();
        j = 0;
        while(((power(g, j * k)*1) != (A[i]*1)) && (j <= fieldSize))
        {
          //cout << "power(g, j * k)*1: " << power(g, j * k)*1 << " - A[i]*1: " << A[i]*1 << endl;
          j++;
        }
        if(j == fieldSize + 1)  return(-1);
        result = result + power(g, j) * power(A.mvar(), i);
        //cout << "power(A[i], k): " << power(A[i], k) << endl;
        //cout << "result: " << result << endl;
      }
      else
      {
        //cout << "not inBaseDomain" << endl;
        result = result + GFPowDown(A[i], k) * power(A.mvar(), i);;
        // power(pol[i], k) * power(A.mvar(), i);
        //pol[i] = GFMap(pol[i], k);
      }
    }
  }
  return(result);
}

static CanonicalForm GFPowUp(const CanonicalForm & A, int k)
// Raises all coefficients of A in the ground field to the power k.
// Used for maping a polynomial to an extension field.
{
  //cout << "A:" << A <<"\n";
  CanonicalForm result = 0;
  int i;
  for(i = 0; i <= degree(A); i++)
  {
    if(!A[i].isZero())
    {
      //cout << "i: " << i << endl;
      //cout << "A[i]: " << A[i] << endl;
      if(A[i].inBaseDomain())
      {
        //cout << "inBaseDomain" << endl;
        result = result + power(A[i], k) * power(A.mvar(), i);
        //cout << "power(A[i], k): " << power(A[i], k) << endl;
        //cout << "result: " << result << endl;
      }
      else
      {
        //cout << "not inBaseDomain" << endl;
        result = result + GFPowUp(A[i], k) * power(A.mvar(), i);;
      }
    }
  }
  //cout << "result:" << result <<"\n";
  return(result);
}

CanonicalForm GFMapUp(const CanonicalForm & A, int k)
// Maps all coefficients of A to the base domain, asumming A is in GF(p^k).
// The current base domain must be GF(p^j), with j a multiple of k.
{
  int p = getCharacteristic();
  int expon = getGFDegree();
  int extExp = expon / k;
    // Assumes that we are using Conway polynomials
  return(GFPowUp(A, extExp));
}

CanonicalForm GFMapDown(const CanonicalForm & A, int k)
// Maps all coefficients of A from the base domain to GF(p^k).
// The current base domain must be GF(p^j), with j a multiple of k.
{
  int p = getCharacteristic();
  int expon = getGFDegree();
  //cout << "Expon: " << expon << endl;
  int extExp = expon / k;
    // Assumes that we are using Conway polynomials
  return(GFPowDown(A, extExp));
}

CanonicalForm newGCD(CanonicalForm A, CanonicalForm B)
// Computes the GCD of two polynomials over a prime field.
// Based on Algorithm 7.2 from "Algorithms for Computer Algebra"
{
  //int debug = 1;
  //if(debug)
  //{
  //  cout << "newGCD input" << endl;
  //  cout << "A: " << A << endl;
  //  cout << "B: " << B << endl;
  //}

  if(A*1 == 0)
  {
    //if(debug)cout << "A Zero!" << endl;
    return(abs(B));
  }
  if(B*1 == 0)
  {
    //if(debug)cout << "B Zero!" << endl;
    return(abs(A));
  }
  if (A.inBaseDomain()) return A.genOne();
  if (B.inBaseDomain()) return B.genOne();

  CanonicalForm pol;
  int p = getCharacteristic();
  int k = 1;
  int fieldSize = p;
  if (CFFactory::gettype() == GaloisFieldDomain)
  {
    k=getGFDegree();
    fieldSize = pow(p, k);
  }

  //if(debug)
  //{
  //  cout << "p: " << p << endl;
  //  cout << "k: " << k << endl;
  //  cout << "fieldSize: " << fieldSize << endl;
  //  cout << "pow: " << pow(p, k) << endl;
  //}

  CFMap M,N;
  compress(A,B,M,N);
  A=M(A);
  B=M(B);
  // Compute the main variable (the largest in A and B)
  Variable x = A.mvar();
  int mv = x.level();
  if(B.mvar().level() > mv)
  {
    x = B.mvar();
    mv = x.level();
  }
  //cout << "main variable " << x << endl << endl;

  // Checks for univariate polynomial
  if ( mv == 1 )
  {
    // This call can be replaced by a faster GCD algorithm
    pol = simpleGCD(A, B); // Univariate GCD
    //cout << "newGCD output 1: " << pol << endl;
    return N(pol);
  }

  CanonicalForm b;
  CFArray bArray(0, fieldSize-1);    // Stores the bs already used
  int sizebArray;

  GFRandom genGF;
  FFRandom genFF;
  int i, j;
  int used;

  CanonicalForm c;    // gcd of the contents
  CanonicalForm C;
  CanonicalForm Cb;    // gcd of Ab and Bb
  CanonicalForm H = 0;      // for Newton Interpolation
  CanonicalForm newtonPoly = 1;  // for Newton Interpolation
  CanonicalForm Cblc;

  // Contents and primparts of A and B
  CanonicalForm contA, contB;    // Contents of A and B
  CanonicalForm Ap, Bp;    // primpart of A and B
  contA = contentWRT(A, mv - 1);
  contB = contentWRT(B, mv - 1);
  c = simpleGCD(contA, contB);    // gcd of univariate polynomials
  Ap = A / contA;
  Bp = B / contB;

  CanonicalForm AL, BL;  // leading coefficients of A and B
  CanonicalForm g;    // gcd of AL and BL
  AL = lcoefWRT(Ap, mv - 1);
  BL = lcoefWRT(Bp, mv - 1);
  g = simpleGCD(AL, BL);    // gcd of univariate polynomials

  CanonicalForm Ab, Bb, gb;  // A, B and g evaluated at b

  int m;          // degree of the gcd of Ab and Bb
  int n = degWRT(Ap, mv - 1);
  if(degWRT(Bp, mv - 1) < n)
  {
    n = degWRT(Bp, mv - 1);
  }

  int fail = 0;
  int goodb;
  // The main loop
  //cout << "start loop\n";
  sizebArray = 0;
  do
  {
    // Searches for a good b.
    // If there are no more possible b, the algorithm fails.
    goodb = 0;
    if(sizebArray >= fieldSize-1)
    {
      // An extension field is needed.
      fail = 1;
    }
    else
    {
      do
      {
        // Searches for a new element of the ground field
        do
        {
          // New element of the ground field
          if(k > 1)
          {
            b = genGF.generate();
          }
          else
          {
            b = genFF.generate();
          }
          //cout << "try(" << sizebArray << "):" << b;
          // Checks if this element was already used
          used = 0;
          for(i = 0; i < sizebArray; i++)
          {
            // Multiplication by 1 is used to prevent a bug which
            // produces b=a^(q-1) as a random element.
            if((bArray[i]*1) == (b*1))
            {
              used = 1;
              //cout << " used\n";
            }
          }
          //if(debug==1)cout << "b: " << b << endl;
        }
        while(used == 1);
        bArray[sizebArray] = b;
        sizebArray++;
        // b must not cancel the gcd of the leading coefficients
        if(g(b, mv) != 0)
        {
          goodb = 1;
          //  cout << " good\n";
        }
        else
        {
          //  cout << " bad";
          if(sizebArray == fieldSize - 1)
          {
            fail = 1;
            //cout << " out of elems " << sizebArray << "tried";
          }
          //cout << "\n";
        }
      }
      while((goodb == 0) && (fail == 0));
    }
    if(fail)
    {
      // Algorithm fails. An extension is needed.

      // Computes the exponent of the ring extension so as to have enough interpolation points.
      int degMax;
      if(totaldegree(A) > totaldegree(B))
      {
        degMax = totaldegree(A);
      }
      else
      {
        degMax = totaldegree(B);
      }
      int expon = 2; // expon <= will not extend the field
      while(pow(fieldSize, expon) < degMax)
      {
        expon++;
      }
      //if(debug)cout << "Not enough elements in the base field. An extension to " << p << "^" << k*expon << " is needed." << endl;
      if(k > 1)
      {
        if(pow(p,k * expon) < (1<<16))
        {
          setCharacteristic(p, k * expon, 'b');
          CanonicalForm P1 = GFMapUp(A, k);
          CanonicalForm P2 = GFMapUp(B, k);
          //cout << "newGCD(mapped):" << P1 << " , " << P2 <<"\n";
          pol = newGCD(P1, P2);
          pol = GFMapDown(pol, k);
          setCharacteristic(p, k, 'a');
          //if(debug)cout << "newGCD output 5: " << pol << endl;
          CanonicalForm temp=N(pol);
          temp/=temp.lc();
          return temp;
        }
        else
        {
          Off(SW_USE_GCD_P);
          CanonicalForm temp=N(gcd( A, B));
          On(SW_USE_GCD_P);
          return temp;
        }
      }
      else
      {
        if(pow(p,k * expon) < (1<<16))
        {
          setCharacteristic(p, k * expon, 'a');
          CanonicalForm P1 = A.mapinto();
          CanonicalForm P2 = B.mapinto();
          pol = newGCD(P1, P2);
          setCharacteristic(p);
          //if(debug)cout << "newGCD output 4: " << pol << endl;
          CanonicalForm temp=N(pol);
          temp/=temp.lc();
          return temp;
        }
        else
        {
          Off(SW_USE_GCD_P);
          CanonicalForm temp=N(gcd( A, B));
          On(SW_USE_GCD_P);
          return temp;
        }
      }
    }
    else
    {
      // Evaluate the polynomials in b
      Ab = Ap(b, mv);
      Bb = Bp(b, mv);
      gb = g(b, mv);
      Cb = newGCD(Ab, Bb);
      //if(debug)cout << "newGCD received: " << Cb << endl;
      m = Cb.degree();

      Cblc = Cb.lc();
      Cb *= gb;
      Cb /= Cblc;
      // Test for unlucky homomorphism
      if(m < n)
      {
        // The degree decreased, we have to start it all over.
        H = Cb;
        newtonPoly = newtonPoly * (x - b);
        n = m;
      }
      else if(m == n)
      {
        // Good b
        H = newtonInterp(b, Cb, newtonPoly, H, x);
        newtonPoly = newtonPoly * (x - b);
        if(lcoefWRT(H, mv - 1) == g)
        {
          C = H / contentWRT(H, mv - 1);  // primpart
          if(fdivides(C, A) && fdivides(C, B))
          {
            //if(debug)cout << "newGCD output 2: " << c * C<< endl;
            return N(c * C);
          }
          else
          {
            if(m == 0)
            {
              //if(debug)cout << "newGCD output 3: " << c << endl;
              return N(c);
            }
          }
        }
      }
    }
  } while(1);
  //cout << "No way to get here!" << endl;
  //cout << "H:" << H << endl;
  // No way to get here!
  return H;
}

#if 0
main()
{
  CanonicalForm A;
  CanonicalForm B;
  CanonicalForm pol;
  Variable x('x'), y('y'), z('z');

  cout << "setCharacteristic(2, 4, 'a')" << endl;
  setCharacteristic(2, 4, 'a');

  int k = getGFDegree();
  int p = getCharacteristic();
  int fieldSize = pow(p, k);
  cout << "p: " << p << endl;
  cout << "GFDegree: " << k << endl;
  cout << "fieldSize: " << fieldSize << endl << endl;

  CanonicalForm g = getGFGenerator();
  //CanonicalForm g = 1;
  A = power((x*y*y-power(x, 7)), 42) * power((power(y, 5) + g*y*y*power(x,13) + g*g), 13);
  B = power((x+y+1), 37) * power((power(y, 5) + g*y*y*power(x,13) + g*g), 13);
  cout << "A: " << A << endl;
  cout << "B: " << B << endl << endl;

  int i;
  CanonicalForm qa;
  CanonicalForm lco;

  for(i = 1; i <= 1; i++){
    pol = newGCD(A * i, B * i);
    lco = pol.lc();
    cout << "new GCD: " << (pol / lco) << endl;
  }

  for(i = 1; i<=1; i++){
    pol = gcd(A*i, B*i);
    lco = pol.lc();
    cout << "old GCD: " << pol / lco << endl;
  }

}
#endif

