/*
   mpn_mulmid-tune.c:  tuning program for integer middle product algorithms

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
   Finds approximate threshold between basecase and karatsuba integer middle
   product algorithms.

   Store this in mpn_smp_kara_thresh, and writes some logging information
   to flog.
*/
void
tune_mpn_smp_kara (FILE* flog, int verbose)
{
   fprintf (flog, "mpn smp kara: ");
   fflush (flog);

   // how long we are willing to wait for each profile run
   const double speed = 0.0001;

   // reset the threshold so karatsuba never accidentally recurses
   ZNP_mpn_smp_kara_thresh = SIZE_MAX;

   size_t thresh;
   const int max_intervals = 200;
   size_t points[max_intervals + 1];
   double score[max_intervals + 1];
   double result[2];
   profile_info_t info;

   // find an upper bound, where karatsuba appears to be safely ahead of
   // basecase
   size_t upper;
   int found = 0;
   for (upper = 2; upper <= 1000 && !found; upper = 2 * upper)
   {
      info->n1 = 2 * upper - 1;
      info->n2 = info->n = upper;

      result[0] = profile (NULL, NULL, profile_mpn_smp_basecase, info, speed);
      result[1] = profile (NULL, NULL, profile_mpn_smp_kara, info, speed);

      if (result[1] < 0.9 * result[0])
         found = 1;
   }

   if (!found)
   {
      // couldn't find a reasonable upper bound
      thresh = SIZE_MAX;
      goto done;
   }

   // subdivide [2, upper] into intervals and sample at each endpoint
   double lower = 2.0;
   double ratio = (double) upper / lower;
   unsigned i;
   for (i = 0; i <= max_intervals; i++)
   {
      points[i] = ceil (lower * pow (ratio, (double) i / max_intervals));
      info->n1 = 2 * points[i] - 1;
      info->n2 = info->n = points[i];
      result[0] = profile (NULL, NULL, profile_mpn_smp_basecase, info, speed);
      result[1] = profile (NULL, NULL, profile_mpn_smp_kara, info, speed);
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
      if (thresh == SIZE_MAX)
         fprintf (flog, "infinity");
      else
         fprintf (flog, "%lu", thresh);
   }
   else
      fprintf (flog, "done");

   fflush (flog);

   ZNP_mpn_smp_kara_thresh = thresh;

   fprintf (flog, "\n");
}



/*
   Finds approximate threshold between karatsuba and fallback integer middle
   product algorithms.

   Store this in mpn_mulmid_fallback_thresh, and writes some logging
   information to flog.
*/
void
tune_mpn_mulmid_fallback (FILE* flog, int verbose)
{
   fprintf (flog, "mpn mulmid fallback: ");
   fflush (flog);

   // how long we are willing to wait for each profile run
   const double speed = 0.0001;

   size_t thresh;
   const int max_intervals = 30;
   size_t points[max_intervals + 1];
   double score[max_intervals + 1];
   double result[2];
   profile_info_t info;

   // find an upper bound, where fallback appears to be safely ahead of
   // karatsuba
   size_t upper = ZNP_mpn_smp_kara_thresh;
   int found = 0;
   for (; upper <= 40000 && !found; upper = 2 * upper)
   {
      info->n1 = 2 * upper - 1;
      info->n2 = info->n = upper;

      result[0] = profile (NULL, NULL, profile_mpn_smp_kara, info, speed);
      result[1] = profile (NULL, NULL, profile_mpn_mulmid_fallback,
                           info, speed);

      if (result[1] < 0.9 * result[0])
         found = 1;
   }

   if (!found)
   {
      // couldn't find a reasonable upper bound
      thresh = SIZE_MAX;
      goto done;
   }

   // subdivide [kara_thresh, upper] into intervals and sample at
   // each endpoint
   double lower = ZNP_mpn_smp_kara_thresh;
   double ratio = (double) upper / lower;
   unsigned i;
   for (i = 0; i <= max_intervals; i++)
   {
      points[i] = ceil (lower * pow (ratio, (double) i / max_intervals));
      info->n1 = 2 * points[i] - 1;
      info->n2 = info->n = points[i];
      result[0] = profile (NULL, NULL, profile_mpn_smp_kara, info, speed);
      result[1] = profile (NULL, NULL, profile_mpn_mulmid_fallback,
                           info, speed);
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
      if (thresh == SIZE_MAX)
         fprintf (flog, "infinity");
      else
         fprintf (flog, "%lu", thresh);
   }
   else
      fprintf (flog, "done");

   fflush (flog);

   ZNP_mpn_mulmid_fallback_thresh = thresh;

   fprintf (flog, "\n");
}


// end of file ****************************************************************
