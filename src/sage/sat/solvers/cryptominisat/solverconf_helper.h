/**
 * @file solverconf_helper.h
 *
 * This file implements a map from option names to C++ class
 * members. This map allows to maintain only one list of options in
 * the whole Cython interface which means that maintainance is easier
 * if CryptoMiniSat's interface should change.
 */

/*****************************************************************************
 *  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  The full text of the GPL is available at:
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include <cryptominisat/Solver/SolverConf.h>

enum sc_type {
  t_int      = 1<<0,
  t_float    = 1<<1,
  t_double   = 1<<2,
  t_Var      = 1<<3,
  t_bool     = 1<<4,
  t_uint32_t = 1<<5,
  t_uint64_t = 1<<6,
};

typedef struct {
  sc_type type;
  const char *name;
  void *target;
  const char *doc;
} sc_entry;

size_t setup_map(sc_entry *entries, CMSat::SolverConf &solver_conf, const size_t n);
