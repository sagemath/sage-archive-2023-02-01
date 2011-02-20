/* IA-32 floating point unit non-inline related functions.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2011 BUGSENG srl (http://bugseng.com)

This file is part of the Parma Polyhedra Library (PPL).

The PPL is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The PPL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02111-1307, USA.

For the most up-to-date information see the Parma Polyhedra Library
site: http://www.cs.unipr.it/ppl/ . */

#include <ppl-config.h>

#if PPL_CAN_CONTROL_FPU && defined(PPL_FPMATH_MAY_USE_SSE) \
  && defined(__i386__) \
  && (defined(__GNUC__) || defined(__INTEL_COMPILER))

#include "fpu.defs.hh"
#include <csetjmp>
#include <csignal>

namespace {

jmp_buf env;

void
illegal_instruction_catcher(int) {
  signal(SIGILL, SIG_DFL);
  longjmp(env, 1);
}

} // namespace

namespace Parma_Polyhedra_Library {

bool have_sse_unit = true;

void
detect_sse_unit() {
  sighandler_t old_handler = SIG_DFL;

  if (setjmp(env)) {
    // We will end up here if sse_get_control() raises SIGILL.
    have_sse_unit = false;
    goto restore_sigill_handler;
  }

  // Install our own signal handler for SIGILL.
  old_handler = signal(SIGILL, illegal_instruction_catcher);
  (void) sse_get_control();
  // sse_get_control() did not raise SIGILL: we have an SSE unit.
  have_sse_unit = true;

 restore_sigill_handler:
  // Restore the previous signal handler for SIGILL.
  signal(SIGILL, old_handler);
}

} // namespace Parma_Polyhedra_Library

#endif // PPL_CAN_CONTROL_FPU && defined(PPL_FPMATH_MAY_USE_SSE) && defined(__i386__) && (defined(__GNUC__) || defined(__INTEL_COMPILER))
