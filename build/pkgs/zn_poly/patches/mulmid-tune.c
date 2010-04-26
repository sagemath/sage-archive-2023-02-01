/*
   mulmid-tune.c:  tuning program for middle product algorithms

   Copyright (C) 2007, 2008, David Harvey

   This file is part of the zn_poly library (version 0.9).

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <math.h>
#include "support.h"
#include "zn_poly_internal.h"
#include "profiler.h"


/*
   For each modulus size, finds approximate threshold between KS4 and fft
   middle product algorithms.

   (Note this needs to be done *after* the KS1/KS2/KS4 multiplication and
   nussbaumer thresholds have been determined, since they are used as
   subroutines.)

   Store these in the global threshold table, and writes some logging
   information to flog.
*/
void
tune_mulmid (FILE* flog, int verbose)
{
   unsigned b;

   fprintf (flog, " KS/FFT mulmid: ");
   fflush (flog);

   // how long we are willing to wait for each profile run
   const double speed = 0.0001;

   // run tuning process for each modulus size
   for (b = 2; b <= ULONG_BITS; b++)
   {
      // thresh for KS4->FFT
      size_t thresh;

      const int max_intervals = 20;
      size_t points[max_intervals + 1];
      double score[max_intervals + 1];

      double result[2];

      profile_info_t info[2];
      info[0]->m = info[1]->m = (1UL << (b - 1)) + 1;

      info[0]->algo = ALGO_MULMID_KS4;
      info[1]->algo = ALGO_MULMID_FFT;

      // find an upper bound, where FFT algorithm appears to be safely
      // ahead of KS4 algorithm
      size_t upper;
      int found = 0;
      for (upper = 45; upper <= 65536 && !found; upper = 2 * upper)
      {
         info[0]->n1 = info[1]->n1 = 2 * upper;
         info[0]->n2 = info[1]->n2 = upper;

         result[0] = profile (NULL, NULL, profile_mulmid, info[0], speed);
         result[1] = profile (NULL, NULL, profile_mulmid, info[1], speed);

         if (result[1] < 0.95 * result[0])
            found = 1;
      }

      if (!found)
      {
         // couldn't find a reasonable upper bound
         thresh = SIZE_MAX;
         goto done;
      }

      // find a lower bound, where KS4 algorithm appears to be safely
      // ahead of FFT algorithm
      size_t lower;
      found = 0;
      for (lower = upper/2; lower >= 2 && !found; lower = lower / 2)
      {
         info[0]->n1 = info[1]->n1 = 2 * lower;
         info[0]->n2 = info[1]->n2 = lower;

         result[0] = profile (NULL, NULL, profile_mulmid, info[0], speed);
         result[1] = profile (NULL, NULL, profile_mulmid, info[1], speed);

         if (result[1] > 1.05 * result[0])
            found = 1;
      }

      if (!found)
      {
         // couldn't find a reasonable lower bound
         thresh = 0;
         goto done;
      }

      // subdivide [lower, upper] into intervals and sample at each endpoint
      double ratio = (double) upper / (double) lower;
      unsigned i;
      for (i = 0; i <= max_intervals; i++)
      {
         points[i] = ceil (lower * pow (ratio, (double) i / max_intervals));
         info[0]->n1 = info[1]->n1 = 2 * points[i];
         info[0]->n2 = info[1]->n2 = points[i];
         result[0] = profile (NULL, NULL, profile_mulmid, info[0], speed);
         result[1] = profile (NULL, NULL, profile_mulmid, info[1], speed);
         score[i] = result[1] / result[0];
      }

      // estimate threshold
      unsigned count = 0;
      for (i = 0; i <= max_intervals; i++)
         if (score[i] > 1.0)
            count++;
      thresh = (size_t)
           ceil (lower * pow (ratio, (double) count / (max_intervals + 1)));

      done:

      if (verbose)
      {
         fprintf (flog, "\nbits = %u, cross to FFT at ", b);
         if (thresh == SIZE_MAX)
            fprintf (flog, "infinity");
         else
            fprintf (flog, "%lu", thresh);
      }
      else
         fprintf (flog, ".");

      fflush (flog);

      tuning_info[b].mulmid_fft_thresh = thresh;
   }

   fprintf (flog, "\n");
}


// end of file ****************************************************************
