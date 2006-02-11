#include <fstream>

#include "power_series.h"


// Construct an empty power series              // *** TO DO: Phase this out!
template <class T>
PowerSeries<T>::PowerSeries() {

  // Make an empty power series
  _precision = -1;
  _default_var = 'q';

}


// Construct a power series of a fixed precision
template <class T>
PowerSeries<T>::PowerSeries(const long precision){

  vector<T> new_series(precision+1, 0);  // Fill the new_series vector with zeros

  _series = new_series;
  _precision = precision;
  _default_var = 'q';

}



// Constructor using "Make(QQ)" to initialize the local conditions
template <class T>
PowerSeries<T>::PowerSeries(const vector<T> & vec){
  _series = vec;
  _precision = _series.size();
  _default_var = 'q';

}


// Empty Destructor
template <class T>
PowerSeries<T>::~PowerSeries(){
  // Do nothing...
}





/////////////////////
// Copy Constructor
/////////////////////
template <class T>
void PowerSeries<T>::operator=(const PowerSeries & source) {

  // Protect against self-assignment
  if (this != &source) {
    _series = source._series;
    _precision = source._precision;
    _default_var = source._default_var;
  }

}







/////////////////////////////////////////////////////////////////////////////
// Return the precision of the power series (so it's O(q^{precision + 1})
///////////////////////////////////////////////////////////////////////////////
template <class T>
long PowerSeries<T>::Precision() const {

  // To Do: We should make a separate empty_flag, and then keep precision as an unsigned long...
  //        but this requires us to filter all routines for the empty power series!

  return _precision;

}



//////////////////////////////////////////////////////
// Return the q^i-th coefficient of the power series
//////////////////////////////////////////////////////
template <class T>
const T & PowerSeries<T>::operator[](const unsigned long i) const {

  // Check we're not out of range...
  if ((i < 0) || (i > _precision)) {
    cout << "Error in PowerSeries<T>::operator[]:  The index " << i << " is out of range! =( " << endl;
    exit(1);
  }

  // Otherwise return the value
  return _series[i];
}

template <class T>
T & PowerSeries<T>::operator[](const unsigned long i) {

  // Check we're not out of range...
  if ((i < 0) || (i > _precision)) {
    cout << "Error in PowerSeries<T>::operator[]:  The index " << i << " is out of range! =( " << endl;
    exit(1);
  }

  // Otherwise return the value
  return _series[i];
}




///////////////////////////////////////////////////////
// Compares when two local_repn_arrays are the same.
////////////////////////////////////////////////////////
template <class T>
bool PowerSeries<T>::operator==(const PowerSeries & B) const {

  // Check they have the same size
  if ((*this)._precision != B._precision)
    return false;

  // Check they have the same entries
  for(long i=0; i < _precision; i++)
    if ((*this)._series[i] != B._series[i])
      return false;

  // If both, then they're equal! =)
  return true;
}



///////////////////////////////////////////////////////
// Compares when two local_repn_arrays are different.
////////////////////////////////////////////////////////
template <class T>
bool PowerSeries<T>::operator!=(const PowerSeries & B) const {

  return !((*this) == B);

}



///////////////////////////////////////////////////////
// Says whether the local condition is empty or not
////////////////////////////////////////////////////////
template <class T>
bool PowerSeries<T>::IsEmpty() const {

  return ((*this)._precision < 0);

}




////////////////////////////////
// Prints the power series conditions
///////////////////////////////
template <class T>
ostream & PowerSeries<T>::Print(ostream & out, const char y) const {

  // Deal with the empty power series
  if ((*this).IsEmpty() == true)
    out << "Empty PowerSeries";
  else {

    // Decide on the variable name
    char x;
    if (y != ' ')
      x = y;
    else
      x = _default_var;


    // Print the constant term
    out << _series[0];

    // Print the first term
    if (_precision >= 1)
      out << " + " << _series[1] << "*" << x;

    // Print all other terms
    if (_precision >= 2)
      for (long i=2; i <= _precision; i++)
	out << " + " << _series[i] << "*" << x << "^" << i;

    // Print the error term
    out << " + O(" << x << "^" << (_precision + 1) << ")";

  }


  return out;
}



////////////////////////////////////////////////////
// This is the overloaded (global) version for <<
////////////////////////////////////////////////////
template <class T>
ostream & operator<<(ostream & out, const PowerSeries<T> & PS) {
  return PS.Print(out);
}



////////////////////
// Define addition
////////////////////
template <class T>
void PowerSeries<T>::operator+(const PowerSeries & B) {

  // Set the precision (dealing with the empty series at the same time)
  _precision = min(_precision, B._precision);

  // Add the new series to the current one
  for(long i=0; i <= _precision; i++)
    _series[i] += B._series[i];

}


//////////////////////
// Define subtrction
//////////////////////
template <class T>
void PowerSeries<T>::operator-(const PowerSeries & B) {

  // Set the precision (dealing with the empty series at the same time)
  _precision = min(_precision, B._precision);

  // Subtract the new series to the current one
  for(long i=0; i <= _precision; i++)
    _series[i] -= B._series[i];

}



//////////////////////////
// Define multiplication
//////////////////////////
template <class T>
void PowerSeries<T>::operator*(const PowerSeries & B) {

  // Set the precision (dealing with the empty series at the same time)
  _precision = min(_precision, B._precision);

  // Make a new vector (of zeros) to hold the multiplication
  vector<T> temp_vec(_precision, 0);

  // Multiply the two series together
  for(long i=0; i <= _precision; i++)
    for(long j=0; j <= i; j++)
      temp_vec[i] += (_series[j] * B._series[i-j]);

  // Replace the current series with the product
  _series = temp_vec;

}




//////////////////////////////////
// Reads the series from a file
//////////////////////////////////
template <class T>
void PowerSeries<T>::ReadSeries(const char* SeriesFilename, const bool use_file_precision) {


  cout << " Reading the file: " << SeriesFilename <<  endl;


  // Try to open the file
  ifstream seriesfile;
  seriesfile.open(SeriesFilename, ios::in);

  // Abort if we fail... =(
  if (! seriesfile.is_open())
    { cout << "ReadSeries Error: Error opening file"; exit(1); }


  // Read the series (assuming the terms are in order of ascending exponents)
  char c;
  T num;
  long pow;
  bool DoneFlag = false;

  pow = -1;

  while ((DoneFlag == false) && (pow < _precision) && (seriesfile.eof() == false)) {

    // cout << "\n Just finished pow = " << pow << endl;

    // Check to see if we're done "O(q^...)"
    seriesfile >> c;  // cout << " Read in char " << c << endl;
    if (c != 'O') {
      seriesfile.putback(c);  // cout << " Put back char " << c << endl;

      // Check to see of there is no leading coefficient (which means the coefficient is 1)
      if (c == 'q') {
	num = 1;
	c = '*';
      }
      else {
	// Get the constant term
	seriesfile >> num;  // cout << " Read in num = " << num << endl;

	// Eat the "+" or "*"
	seriesfile >> c;  // cout << " Read in char " << c << endl;
      }


      if (c == '+')
	_series[0] = num;
      else {

	if (c != '*') {
	  cout << " ReadSeries Error1: Unexpected input character " << c << " found." << endl;
	  exit(1);
	}


	/*
	// DIAGNOSTIC
	cout << " Checking for a linear term " << endl;
	*/


	// Check if there is a linear term

	// Eat the "q +" or "q^"
	seriesfile >> c;  // cout << " Read in char " << c << endl;
	seriesfile >> c;  // cout << " Read in char" << c << endl;

	if (c == '+')
	  _series[1] = num;
	else {

	  if (c != '^') {
	    cout << " ReadSeries Error2: Unexpected input character " << c << " found." << endl;
	    exit(1);
	  }


	  /*
	  // DIAGNOSTIC
	  cout << " Checking for higher terms " << endl;
	  */


	  // See where to put the term if it has an exponent
	  seriesfile >> pow;  // cout << " Read in pow = " << pow << endl;


	  // Eat the final "+"
	  seriesfile >> c;  // cout << " Read in char (after the power) " << c << endl;


	  // Store the coefficient if it's in range
	  if ((pow > 0) && (pow <= _precision))
	    _series[pow] = num;
	  else {
	    /*
	    cout << " ReadSeries Warning: Exponent " << pow << " of series is out of range (0 ... " << _precision << ")." << endl;
	    cout << " Disregarding exponent " << pow << "." << endl;
	    */
	  }
	}
      }

    }
    else {

      // Read the precision + 1
      unsigned long file_precision;
      seriesfile >> c;  // Reads '('
      seriesfile >> c;  // Reads 'q'
      seriesfile >> c;  // Reads ')' or '^'
      if (c == ')')
	file_precision = 1;
      else {
	seriesfile >> file_precision;
	file_precision--;
      }

      // If use_file_precision == true, then set this to be the new precision,
      // otherwise check it's compatible with the expected precision
      if (use_file_precision == true)
	_precision = file_precision;
      else
	if (_precision > file_precision) {
	  cout << "Error in ReadFile():  The file precision is not as large as we expected! =( " << endl;
	  exit(1);
	}

      // In any case, we're done reading the file! =)
      DoneFlag = true;
    }

  }


  // Close the file
  seriesfile.close();

}








////////////////////////////////
// Writes the series to a file
////////////////////////////////
template <class T>
void PowerSeries<T>::WriteSeries(const char* SeriesFilename, const char y) {

  // Try to open the file
  ofstream seriesfile;
  seriesfile.open(SeriesFilename, ios::out);

  // Abort if we fail... =(
  if (! seriesfile.is_open()) {
    cout << "WriteSeries Error: Error opening file" << endl;
    cout << "  Tried to open filename: " << SeriesFilename << endl << endl;
    exit (1);
  }

  // Write the file
  (*this).Print(seriesfile, y);

  // Close the file
  seriesfile.close();

}












