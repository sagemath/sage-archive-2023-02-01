
///////////////////////////////////////////////
// Finds the minimum of two rational numbers //
///////////////////////////////////////////////

mpq_class Minimum(mpq_class a, mpq_class b) {

  if (a > b)
    return b;
  else
    return a;
}


/////////////////////////////////////////////////////////
// Finds the minimum of a valarray of rational numbers //
/////////////////////////////////////////////////////////

mpq_class Minimum(valarray <mpq_class> v) {

  mpq_class min_elt;
  min_elt = v[0];    // This will give an error if v is empty... which is good. =)

  for (size_t i = 1; i < v.size(); i++)
    if (min_elt > v[i])
      min_elt = v[i];

  return min_elt;
}
