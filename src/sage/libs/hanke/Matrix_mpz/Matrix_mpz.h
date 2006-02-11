#ifndef MATRIX_MPZ_H
#define MATRIX_MPZ_H

#include "Matrix_mpz_header.h"

#include "../LCM.h"   // Do we need this??
#include "../Utilities/string_utils.h"    // Needed for the filenames in Eisenstein_Bound.cc

#include "output.h"   // Do we need this??


// Basic Matrix Stuff
#include "Matrix_mpz.cc"


// Some Extras we need
#include "../GMP_class_extras/mpz_class_extras.h"
#include "../GMP_class_extras/vectors.h"
#include "../GMP_class_extras/mpq_stuff.h"
#include "../GMP_class_extras/Bernoulli.h"

/* NOT INCLUDED IN SAGE
#include "../Siegel_Diagnostic/siegel_diagnostic.h"   // Needed for the ReadSeries() routine used in "Eisenstein_Bound.cc"
*/


// Local Stuff
#include "Local_Normal.cc"
#include "Matrix_mpz_Extras.cc"

#include "Local_Density_Front.cc"
#include "CountLocal.cc"
#include "CountLocal2.cc"
#include "Local_Density_Congruence.cc"

#include "Local_Invariants.cc"
#include "Local_Constants.cc"


/* NOT INCLUDED IN SAGE
#include "Siegel_Product.cc"
#include "Local_Diagnostic.cc"
#include "Eisenstein_Bound.cc"
*/

#endif
