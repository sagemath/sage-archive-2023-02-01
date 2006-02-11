#include "string_utils_header.h"

// Converts a long to a string
string MakeString(const long & num) {
  ostringstream new_string;
  new_string << num;
  return new_string.str();
}


// Converts an unsigned long to a string
string MakeString(const unsigned long & num) {
  ostringstream new_string;
  new_string << num;
  return new_string.str();
}


// Converts a double to a string
string MakeString(const double & num) {
  ostringstream new_string;
  new_string << num;
  return new_string.str();
}


// Converts an mpz_class to a string
string MakeString(const mpz_class & num) {
  ostringstream new_string;
  new_string << num;
  return new_string.str();
}


// Converts an mpq_class to a string
string MakeString(const mpq_class & num) {
  ostringstream new_string;
  new_string << num;
  return new_string.str();
}


// -------------------------------------------------------------


// Makes a string (for use in a temporary filename) based on the current host and PID number
string MakeUniqueFilename() {

  // Get the current hostname
  char host[100];
  gethostname(host, 100);

  ostringstream new_stringstream;
  new_stringstream << host << "__pid=" << getpid();

  return new_stringstream.str();

}
