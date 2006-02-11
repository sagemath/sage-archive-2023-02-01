

// Print the vector V
void PrintV(vector<long> V) {
  cout << "[ ";
  for (size_t i = 0; i < V.size() - 1; i++)
    cout << V[i] << ", ";
  cout << V[V.size()-1] << " ]" << endl;
}



// Print the first n entries of a vector V
void PrintHeadV(vector<long> V, unsigned long n) {
  cout << "[ ";
  unsigned long end = n;
  if (V.size() < n)
    end = V.size();

  for (size_t i = 0; i < end; i++)
    cout << V[i] << ", ";
  cout << " ..." << endl;
}



// Print the first n entries of a vector V
void PrintHeadV(vector<double> V, unsigned long n) {
  cout << "[ ";
  unsigned long end = n;
  if (V.size() < n)
    end = V.size();

  for (size_t i = 0; i < end; i++)
    cout << V[i] << ", ";
  cout << " ..." << endl;
}



// Print the last n entries of a vector V
void PrintTailV(vector<long> V, unsigned long n) {
  cout << "... ";
  unsigned long start = 0;
  if (V.size() > n)
    start = V.size() - n;

  for (size_t i = start; i < V.size() - 1; i++)
    cout << V[i] << ", ";
  cout << V[V.size()-1] << " ]" << endl;
}


// Print the last n entries of a vector V
void PrintTailV(vector<double> V, unsigned long n) {
  cout << "... ";
  unsigned long start = 0;
  if (V.size() > n)
    start = V.size() - n;

  for (size_t i = start; i < V.size() - 1; i++)
    cout << V[i] << ", ";
  cout << V[V.size()-1] << " ]" << endl;
}



// Reads in a vector of longs
vector<long> ReadVector_long(const char* filename) {

    vector<long> Vec;
    long num;
    char c;

    // Open the file
    ifstream infile;
    infile.open(filename, ios::in);

    // Check that the file opened correctly
    if (! infile.is_open())
      { cout << "Error opening file"; exit (1); }


    // Read in the primes
    do {
      infile >> num;
      if (infile.fail() == false)
	Vec.push_back(num);
      infile >> c;   // Discard the comma or find the end of the file...
      //      cout << " Read the number " << num << endl;
    }  while (infile.eof() == false);


    infile.close();

    return Vec;
  }


// Reads in a vector of doubles
vector<double> ReadVector_double(const char* filename) {

    vector<double> Vec;
    double num;
    char c;

    // Open the file
    ifstream infile;
    infile.open(filename, ios::in);

    // Check that the file opened correctly
    if (! infile.is_open())
      { cout << "Error opening file"; exit (1); }


    // Read in the primes
    do {
      infile >> num;
      if (infile.fail() == false)
	Vec.push_back(num);
      infile >> c;   // Discard the comma or find the end of the file...
      //      cout << " Read the number " << num << endl;
    }  while (infile.eof() == false);

    infile.close();

    return Vec;
  }



