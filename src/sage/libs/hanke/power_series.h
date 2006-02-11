#if !defined(POWER_SERIES_H)
#define POWER_SERIES_H
using namespace std;
#include <vector>
#include <iostream>


template <class T>
class PowerSeries {

public:
  PowerSeries();  // Constructor
  PowerSeries(const long precision);
  PowerSeries(const vector<T> & vec);  // Constructor
  ~PowerSeries();  // Destructor

  void operator=(const PowerSeries & source);  // Copy Constructor


  long Precision() const;
  const T & operator[](const unsigned long i) const;
  T & operator[](const unsigned long i);


  bool operator==(const PowerSeries & B) const;  // Compares when two local_repn_arrays are the same.
  bool operator!=(const PowerSeries & B) const;  // Compares when two local_repn_arrays are the same.
  bool IsEmpty() const;


  ostream & Print(ostream & out, const char y = ' ') const;   // Prints the power series


  void operator+(const PowerSeries & B);
  void operator-(const PowerSeries & B);
  void operator*(const PowerSeries & B);


  void ReadSeries(const char* SeriesFilename, const bool use_file_precision = false); // Reads the series from a file
  void WriteSeries(const char* SeriesFilename, const char y = ' ');                    // Writes the series to a file



private:
  vector<T> _series;
  long _precision;     // Warning: This must be signed so we know when it's empty... (_precision == -1)
  char _default_var;

};




// Here are some extras:
// ---------------------
template <class T>
ostream & operator<<(ostream & out, const PowerSeries<T> & PS);



#include "power_series.cc"

// -------------------------------------------------------------------------------------

#endif






