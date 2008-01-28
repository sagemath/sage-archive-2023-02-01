/*
  Copyright 2005, 2006, 2007 David Cadé, Damien Stehlé.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

  This program implements ideas from the paper "Floating-point LLL Revisited",
  by Phong Nguyen and Damien Stehlé, in the Proceedings of Eurocrypt'2005,
  Springer-Verlag; and was partly inspired by Shoup's NTL library:
  http://www.shoup.net/ntl/

*/

#ifndef FPLLL_H
#define FPLLL_H

#if defined(__sun)
#include <ieeefp.h>
extern "C" long double ldexpl(long double x, int exp);
#define NAN __builtin_nanf("")
#endif

#include <mpfr.h>
#include <gmp.h>
#include "dpe.h"

#include "nr.h"
#include "matrix.h"
#include "proved.h"
#include "heuristic.h"
#include "fast.h"
#include "util.h"
#include "wrapper.h"
#include "fast_earlyred.h"
#include "heuristic_early_red.h"

#endif
