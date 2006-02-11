

#if !defined(OUTPUT_H)
#define OUTPUT_H


// Some silly vector/set output operators:
// ---------------------------------------

// Define the << operator for the vector<T> type
template <class T>
ostream & operator<<(ostream & out, const vector<T> & v) {
  out << "[ ";
  for(long i=0; i < v.size(); i++)
    out << v[i] << " ";
  out << " ]";
  return out;
}




// Define the << operator for the set<T> type
template <class T>
ostream & operator<<(ostream & out, const set<T> & s) {

  // Deal with the empty set
  if (s.empty() == true) {
    out << "{ }";
    return out;
  }


  // Deal with non-empty sets
  out << "{ ";
  typename set<T>::iterator i, j;
  //  set<T>::iterator i;                Note: The compiler complains about this.  See Garrett's 1/11/05 e-mail. =)
  for(i = s.begin(); i != s.end(); i++) {
    out << *i;
    j=i;  j++;
    if (j == s.end())
      out << " }";
    else
      out << ", ";
  }

  return out;
}


// Define the >> operator for the set<T> type
template <class T>
istream & operator>>(istream & in, set<T> & s) {

  // Local variables
  set<T> new_set;
  T elt;
  char c;


  // Read the opening brace '{'
  in >> c;

  // Check for the empty set
  in >> c;
  if (c != '}') {

    // If not, then replace the last character on the stream
    in.putback(c);

    // Read in the elements and discard the ',' or '}'.
    do {
      in >> elt;
      if (in.fail() == false)
	new_set.insert(elt);
      in >> c;   // Discard the comma or find the end of the file...
      //      cout << " Read c = " << c << endl;
      //      cout << " Read the number " << num << endl;
    }  while ((c != '}') && (in.eof() == false));

    // Error if we hit the end of the stream!!!
    if (in.eof() == true) {
      cout << " Error in operator >> for sets: Hit the end of the file before reading the closing brace '}' ! =( " << endl;
      assert(0==1);
    }

  }


  // Return the set and the stream
  s = new_set;
  return in;

}









// ============================================================================================


// Include the FileInput() & FileOutput() strict I/O routines
#include "../strict_IO.h"


#endif



