/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Modified copy operators -- defined for each class since we're not allowed to override the template! //
/////////////////////////////////////////////////////////////////////////////////////////////////////////


template<>
valarray<size_t> & valarray<size_t>::operator=(const valarray<size_t> & v) {

  // Protect against self-assignment
  if (this != &v) {
    (*this).resize(v.size());
    for(size_t i=0; i<v.size(); i++)
      (*this)[i] = v[i];
  }

  return *this;
}


template<>
valarray<mpz_class> & valarray<mpz_class>::operator=(const valarray<mpz_class> & v) {

  // Protect against self-assignment
  if (this != &v) {
    (*this).resize(v.size());
    for(size_t i=0; i<v.size(); i++)
      (*this)[i] = v[i];
  }

  return *this;
}

template<>
valarray<mpq_class> & valarray<mpq_class>::operator=(const valarray<mpq_class> & v) {

  // Protect against self-assignment
  if (this != &v) {
    (*this).resize(v.size());
    for(size_t i=0; i<v.size(); i++)
      (*this)[i] = v[i];
  }

  return *this;
}

mpq_class Minimum(const valarray<mpq_class> & vec) {

  // Protect against the empty valarray
  if (vec.size() == 0) {
    cout << "Error in Minimum: The valarray in empty! =(" << endl;
    exit(1);
  }


  // Find the minimum of all of the entries
  mpq_class m;
  m = vec[0];
  for(unsigned long i=1; i < vec.size(); i++)
    m = min(vec[i], m);


  // Return the minimum
  return m;

}
