



//////////////////////
// File I/O for long
//////////////////////
void FileOutput(const long & elt, ostream & out)
{ out << elt; }

void FileInput(long & elt, istream & in)
{ in >> elt; }


///////////////////////////////
// File I/O for unsigned long
///////////////////////////////
void FileOutput(const unsigned long & elt, ostream & out)
{ out << elt; }

void FileInput(unsigned long & elt, istream & in)
{ in >> elt; }


////////////////////////
// File I/O for double
////////////////////////
void FileOutput(const double & elt, ostream & out)
{ out << elt; }

void FileInput(double & elt, istream & in)
{ in >> elt; }


///////////////////////////
// File I/O for mpz_class
///////////////////////////
void FileOutput(const mpz_class & elt, ostream & out)
{ out << elt; }

void FileInput(mpz_class & elt, istream & in)
{ in >> elt; }


///////////////////////////
// File I/O for mpq_class
///////////////////////////
void FileOutput(const mpq_class & elt, ostream & out)
{ out << elt; }

void FileInput(mpq_class & elt, istream & in)
{ in >> elt; }




// ===================================================================================



// File I/O Routines for vector<T>, and some special cases:
// --------------------------------------------------------


//////////////////////////////////////
// Writes a vector<T> to the ostream
//////////////////////////////////////
template <class T>
void FileOutput(const vector<T> & vec, ostream & out) {

  // Run through the vector
  out << " [ ";
  for(size_t i = 0; i < vec.size(); i++) {
    FileOutput(vec[i], out);
    if (i+1 < vec.size())
      out << " , ";
  }
  out << " ] ";

}


///////////////////////////////////////
// Reads a vector<T> from the istream -- replacing the previous vector<T>
///////////////////////////////////////
template <class T>
void FileInput(vector<T> & vec, istream & in) {

  // Read the entries into a temporary vector first
  vector<T> tmp_vec;
  char ch;
  T element;


  // Read the opening " [ "
  in >> ch;
  assert( ch == '[');

  // Quick check that the matrix isn't empty
  // (if so, we'll never enter the loop!)
  in >> ch;
  in.putback(ch);

  // Search for the closing bracket
  while (ch != ']') {

    // Read the next entry and append it to the vector
    FileInput(element, in);
    tmp_vec.push_back(element);

    // Read the separator (',' or ']')
    in >> ch;
    assert((ch == ',') || (ch == ']'));

    // Put back the ']', since it's read at the end
    if (ch == ']')
      in.putback(ch);

  }

  // Read the closing " ] "
  in >> ch;
  assert( ch == ']');

  /*
  // DIAGNSOTIC
  cout << " Started with the vector " << vec << endl;
  cout << " Read the vector " << tmp_vec << endl;
  */

  // Copy the temporary vector for output
  vec = tmp_vec;

  /*
  // DIAGNSOTIC
  cout << " Returned the vector " << vec << endl;
  */




}


