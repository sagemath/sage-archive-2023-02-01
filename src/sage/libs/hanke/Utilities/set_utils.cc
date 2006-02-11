


// Reads in a set of mpz_class
set<mpz_class> ReadSet_mpz_class(const char* filename) {

    set<mpz_class> new_set;
    mpz_class num;
    char c;

    // Open the file
    ifstream infile;
    infile.open(filename, ios::in);

    // Check that the file opened correctly
    if (! infile.is_open())
      { cout << "Error opening file"; exit (1); }


    // Read the opening brace '{'
    infile >> c;


    // Read in the elements and discard the ',' or '}'.
    do {
      infile >> num;
      if (infile.fail() == false)
	new_set.insert(num);
      infile >> c;   // Discard the comma or find the end of the file...
      //      cout << " Read the number " << num << endl;
    }  while (infile.eof() == false);


    // Close the file
    infile.close();

    return new_set;
  }

