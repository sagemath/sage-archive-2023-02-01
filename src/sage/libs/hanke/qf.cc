#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <gmp.h>
#include <fstream>
#include <math.h>
#include <stdlib.h>   // For ultoa() in Matrix_mpz::_GetEisData()
#include <valarray>
#include <time.h>
#include <string.h>

// Additional Libraries from the Older Project...
#include <list>
#include <valarray>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

// Put the source includes here..
// #include "290_Code/GMP_class_extras/mpz_class_extras.h"

#include "Matrix_mpz/Matrix_mpz.h"

// #include "290_Code/Utilities/string_utils.h"
// #include "290_Code/Utilities/file_utilities.h"

// Special code for Python interface.

mpz_class string_to_mpz(const char* s) {
  return mpz_class(s);
}

string mpz_to_string(const mpz_class& x) {
  ostringstream tmp;
  tmp << x;
  return tmp.str();
}

string mpz_to_hex(const mpz_class& x) {
  char* s;
  s = mpz_get_str(NULL, 16, x.get_mpz_t());
  string t(s);
  free(s);
  return t;
}

mpq_class string_to_mpq(const char* s) {
  return mpq_class(s);
}

string mpq_to_string(const mpq_class& x) {
  ostringstream tmp;
  tmp << x;
  return tmp.str();
}

