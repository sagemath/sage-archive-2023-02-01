/* sys/infnan.c
 *
 * Copyright (C) 2001, 2004, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <math.h>

#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_nan (void)
{
  return gsl_fdiv (0.0, 0.0);
}

double gsl_posinf (void)
{
  return gsl_fdiv (+1.0, 0.0);
}

double gsl_neginf (void)
{
  return gsl_fdiv (-1.0, 0.0);
}


int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

#if defined(_MSC_VER) /* Microsoft Visual C++ */
#include <float.h>
int
gsl_isnan (const double x)
{
  return _isnan(x);
}

int
gsl_isinf (const double x)
{
  int fpc = _fpclass(x);

  if (fpc == _FPCLASS_PINF)
    return +1;
  else if (fpc == _FPCLASS_NINF)
    return -1;
  else
    return 0;
}

int
gsl_finite (const double x)
{
  return _finite(x);
}
#else

# if HAVE_DECL_ISFINITE
int
gsl_finite (const double x)
{
  return isfinite(x);
}
# elif HAVE_DECL_FINITE
int
gsl_finite (const double x)
{
  return finite(x);
}
# elif HAVE_IEEE_COMPARISONS
int
gsl_finite (const double x)
{
  const double y = x - x;
  int status = (y == y);
  return status;
}
# endif

# if HAVE_DECL_ISNAN
int
gsl_isnan (const double x)
{
  return isnan(x);
}
#elif HAVE_IEEE_COMPARISONS
int
gsl_isnan (const double x)
{
  int status = (x != x);
  return status;
}
# endif

int
gsl_isinf (const double x)
{
  if (! gsl_finite(x) && ! gsl_isnan(x))
    {
      return (x > 0 ? +1 : -1);
    }
  else
    {
      return 0;
    }
}

#endif

