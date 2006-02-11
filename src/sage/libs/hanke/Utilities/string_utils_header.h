#ifndef STRING_UTILS_HEADER_H
#define STRING_UTILS_HEADER_H

using namespace std;
#include <sstream>
#include <string>
#include <gmp.h>
#include <gmpxx.h>

// Converts a long to a string
string MakeString(const long & num);

// Converts an unsigned long to a string
string MakeString(const unsigned long & num);

// Converts a double to a string
string MakeString(const double & num);

// Converts an mpz_class to a string
string MakeString(const mpz_class & num);

// Converts an mpq_class to a string
string MakeString(const mpq_class & num);

// ----------------------------------------

// Makes a string (for use in a temporary filename) based on the current host and PID number
string MakeUniqueFilename();

#endif
