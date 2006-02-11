from misc import (alarm, srange, xsrange,
                  cputime, verbose, set_verbose, set_verbose_files,
                  get_verbose_files, unset_verbose_files, get_verbose,
                  version, add, union, uniq, powerset, exists, forall,
                  random_sublist, mul, prod, walltime, generic_cmp,
                  repr_lincomb, tmp_dir,
                  DOT_SAGE, SAGE_ROOT, SAGE_URL, SAGE_DB, SAGE_TMP,
                  is_32_bit, is_64_bit, upgrade)

from defaults import set_default_variable_name

from preparser import preparse

from sage_eval import sage_eval

from pyrex import pyrex

from persist import save, load, dumps, loads, db, db_save

from functional import *

from latex import latex, view

from constants import (pi, e, NaN, golden_ratio, log2, euler_gamma, catalan,
                       khinchin, ConstantRing)

from trace import *

##########################################################################
def benchmark(n=-1):
    """
    Run a well-chosen range of SAGE commands and record the time it
    takes for each to run.

    INPUT:
        n -- int (default: -1) the benchmark number; the default
             of -1 runs all the benchmarks.
    OUTPUT:
        list -- summary of timings for each benchmark.
    """
    import sage.misc.benchmark
    return sage.misc.benchmark.benchmark(n)


