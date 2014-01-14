#include "functions.h"

// ------- auxiliary functions ----------------------------
int factorial(int n)
{
  int result=1;
  for (int i=1; i<=n; i++)
    result*=i;
  return(result);
}


// binomial works only well for "normal" values
int binomial(int n, int D)
{
  int d=D;
  if (d>n/2) d=n-d;
  int result=1;
  for (int i=0; i<d; i++)
    result*=n-i;
  for (int i=1; i<=d; i++)
    result/=i;
  return(result);
}
