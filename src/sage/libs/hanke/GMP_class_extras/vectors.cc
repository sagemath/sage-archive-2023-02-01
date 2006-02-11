//#include "vec_mat.h"
#include "vectors.h"
#include "valarray_extras.h"
#include "mpz_class_extras.h"



// Notes:
//    - Sets will mean vectors assumed to be in increasing order, with
//        no duplicates (i.e. strictly increasing).
//    - Routines ending in the word "Ordered" assume that the input is an
//        ordered (i.e. strictly increasing) vector, also referred to as a set.
//        The entries should also be positive integers, but this may not be
//        enforced (I need to look).
//
//    - Types of returns:
//      -----------------
//        Tuples:
//            VectorAppend - append as tuples
//            MakeVector - makes a tuple
//            VectorRemove - makes a tuple
//        Sets:
//            OrderVector - makes a set from a tuple
//            VectorUnion - makes a set
//            VectorIntersection - makes a set
//            VectorComplement - makes a set
//        Boolean (long):
//            IsDisjointOrdered - 0 or 1
//            CheckVectorIntersection - 0 or 1
//            CheckVectorOrdered - 0 or 1




//////////////////////////////////////////////////////////////////////
// Checks if two ordered (strictly increasing) vectors are disjoint //
//////////////////////////////////////////////////////////////////////

bool IsDisjointOrdered(const valarray<int> & v, const valarray<int> & w)
{
  // If either is empty, they're disjoint
  if ((v.size() == 0) || (w.size() == 0))
    return true;

  // Look for a matching entry
  int i, j;
  i = 0; j = 0;  // Start at beginning of v and w
  while ((i < v.size()) && (j < w.size()))
    if (w[j] < v[i])
      j++;
    else if (v[i] < w[j])
      i++;
    else
      return false;  // If they're equal, then the vectors aren't disjoint.

  // If there is no match, then the vectors are disjoint.
  return true;
}


/////////////////////////////////////////////////
// Gives a vector by an arithmetic progression //
/////////////////////////////////////////////////

valarray<int> MakeVector(int length, int start, int increment)
{
  valarray<int> v(length);

  for(int i=0; i<length; i++)
    v[i] = start + i * increment;

  return v;
}


//////////////////////////////////////////////////////////////////////////////////
// Reorders the vector v in strictly increasing order (so no repetitions occur) //
//////////////////////////////////////////////////////////////////////////////////
valarray<int> OrderVector(const valarray<int> & v)
{
  int n, i, j, temp;
  bool done;
  n = v.size();

  // Copy v into v1
  valarray<int> v1;
  v1.resize(n);
  for(long i=0; i<n; i++)
    v1[i] = v[i];


  // Check if v is empty
  if (n == 0)
    return v1;

  // Swap neighboring entries until it's ordered...
  do {
    done = true;
    for (i=0; i < n-1; i++)
      if (v1[i] > v1[i+1]) {
	temp = v1[i];
	v1[i] = v1[i+1];
	v1[i+1] = temp;
	done = false;
      }
  } while (done == 0);


  valarray<int> w, final;
  w.resize(n);

  // Remove duplicate entries
  w[0] = v[0];
  i=1; j=0;
  while (i < n) {
    if (v1[i] > w[j]) {
      j++;
      w[j] = v1[i];
    }
    i++;
  }

  // Resize the final valarray
  final.resize(j+1);
  for (i=0; i <= j; i++)
    final[i] = w[i];

  return(final);
}



//////////////////////////////////////////////
// Creates the union (v + w) of two vectors //
//////////////////////////////////////////////
valarray<int> VectorAppend(const valarray<int> & v, const valarray<int> & w)
{
  int m, n;
  m = v.size();
  n = w.size();

  valarray<int> u;
  u.resize(m+n);

  int i;
  for(i=0; i<m; i++)
    u[i] = v[i];
  for(i=0; i<n; i++)
    u[m+i] = w[i];

  return(u);
}


//////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) union of two vectors //
//////////////////////////////////////////////////////////////////
valarray<int> VectorUnion(const valarray<int> & v, const valarray<int> & w)
{
  return(OrderVector(VectorAppend(v,w)));
}


/////////////////////////////////////////////////////////////////////////////////////
// Finds their ordered union, assuming they are both ordered (strictly increasing) //
/////////////////////////////////////////////////////////////////////////////////////
valarray<int> VectorUnionOrdered(const valarray<int> & v, const valarray<int> & w)
{
  // This could be done more quickly using the ordering,
  // but since out vectors are small, we'll leave this for now... =)
  return(OrderVector(VectorAppend(v,w)));
}


/////////////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) intersection of two vectors //
/////////////////////////////////////////////////////////////////////////
valarray<int> VectorIntersection(const valarray<int> & v, const valarray<int> & w)
{
  int m, n;
  m = v.size();
  n = w.size();

  valarray<int> u;
  u.resize(m+n);

  int i,j,k;
  k=0;
  for(i=0; i<m; i++)
    for(j=0; j<n; j++)
      if (v[i] == w[j]) {
	u[k] = v[i];
	k++;
      }

  // Now trim it to the right size
  valarray<int> final;
  final.resize(k);
  for(i=0; i<k; i++)
    final[i] = u[i];

  return(OrderVector(final));
}


/////////////////////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) intersection of two ordered vectors //
/////////////////////////////////////////////////////////////////////////////////
valarray<int> VectorIntersectionOrdered(const valarray<int> & v, const valarray<int> & w)
{
  // This could be done more quickly using the ordering,
  // but since out vectors are small, we'll leave this for now... =)
  return(VectorIntersection(v,w));
}


////////////////////////////////////
// Remove the entries in w from v //
////////////////////////////////////
valarray<int> VectorRemove(const valarray<int> & v, const valarray<int> & w)
{
  // Deal with the empty vector cases first
  if ((v.size() == 0) || (w.size() == 0))
    return v;

  // Otherwise remove the necessary components
  int m, n;
  m = v.size();
  n = w.size();

  valarray<int> u;
  u.resize(m);

  int i,j, ind;
  ind = 0;
  bool flag;
  for(i=0; i<m; i++) {
    flag = false;

    // Check if v(i) appears in w
    for(j=0; j<n; j++) {
      if (v[i] == w[j])
 	flag = true;
    }

    // Write the entry if it doesn't appear in w
    if (flag == false)
      {
	u[ind] = v[i];
	ind++;
      }
  }

  // Resize the answer
  valarray<int> final;
  final.resize(ind);
  for(i=0; i<ind; i++)
    final[i] = u[i];

  return(final);
}


//////////////////////////////////////////////////////////////
// Returns an ordered vector with the entries of v not in w //
//////////////////////////////////////////////////////////////
valarray<int> VectorComplement(const valarray<int> & v, const valarray<int> & w)
{
  return(OrderVector(VectorRemove(v,w)));
}


///////////////////////////////////////////////////////////////
// Returns an ordered vector with the entries of v not in w, //
// assuming both v and w are ordered (strictly increasing)   //
///////////////////////////////////////////////////////////////
valarray<int> VectorComplementOrdered(const valarray<int> & v, const valarray<int> & w)
{
  // This could be done more quickly using the ordering,
  // but since our vectors are small, we'll leave this for now... =)
  return(OrderVector(VectorRemove(v,w)));
}


/////////////////////////////////////////
// Checks if the two vectors intersect //
/////////////////////////////////////////
bool CheckVectorIntersection(const valarray<int> & v, const valarray<int> & w)
{
  if (VectorUnion(v,w).size() == v.size() + w.size())
    return false;  // No intersection
  else
    return true;  // Some intersection
}


///////////////////////////////////////////////////
// Checks if the two (ordered) vectors intersect //
///////////////////////////////////////////////////
bool CheckVectorIntersectionOrdered(const valarray<int> & v, const valarray<int> & w)
{
  if (VectorUnionOrdered(v,w).size() == v.size() + w.size())
    return false;  // No intersection
  else
    return true;  // Some intersection
}


/////////////////////////////////////////////////
// Checks if v is a strictly increasing vector //
/////////////////////////////////////////////////
bool CheckVectorOrdered(const valarray<int> & v)
{
  int n, i;
  n = v.size();

  for(i=0; i<n-1; i++)
    if (v[i] >= v[i+1])
      return false;               // Return 0 if out of order

  return true;                   // otherwise return 1
}




//////////////////////////////////////////
// Finds the GCD of a vector of numbers //
//////////////////////////////////////////

mpz_class GCD(const valarray<mpz_class> & v)
{
  long n, i;
  n = v.size();

  if (n==1)
    return v[0];

  if (n>1) {
    mpz_class current_GCD;
    current_GCD = v[0];
    for (i=0; i<n-1; i++)
      current_GCD = GCD(current_GCD, v[i+1]);
    return current_GCD;
  }

  cerr << "\n Error in Vector GCD routine: The vector has length n = "
       << n << endl;
  abort();
}




