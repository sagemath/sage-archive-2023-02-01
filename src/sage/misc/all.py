from misc import (alarm, srange, xsrange, sxrange, getitem,
                  cputime, verbose, set_verbose, set_verbose_files,
                  get_verbose_files, unset_verbose_files, get_verbose,
                  version, banner, add, union, uniq, powerset, exists, forall,
                  random_sublist, mul, prod, walltime, generic_cmp,
                  repr_lincomb, tmp_dir, tmp_filename,
                  DOT_SAGE, SAGE_ROOT, SAGE_URL, SAGE_DB, SAGE_TMP,
                  is_32_bit, is_64_bit, newton_method_sizes)

from attach import attach

from profiler import Profiler

from mrange import xmrange, mrange

from bug import bug

from dist import install_scripts

from darcs import darcs_src, darcs_doc, darcs_scripts

from hg import hg_sage, hg_doc, hg_scripts, hg_extcode

from package import install_package, optional_packages, upgrade

from pager import pager

from sagedoc import search_sage

from getusage import top, get_memory_usage

from log import log_html, log_dvi, log_html_mathml

from defaults import set_default_variable_name

from preparser import preparse

from sage_eval import sage_eval

from pyrex import pyrex

from persist import save, load, dumps, loads, db, db_save

from func_persist import func_persist

from functional import *

from latex import latex, view, lprint

from mathml import mathml

from trace import *

from cachefunc import CachedFunction

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


class logstr(str):
    def __repr__(self):
        return self

    def _latex_(self):
        #return "\\begin{verbatim}%s\\end{verbatim}"%self
        if not '#' in self:
         delim = '#'
        elif not '@' in self:
         delim = '@'
        elif not '~' in self:
         delim = '~'
        return r"""\verb%s%s%s"""%(delim, self.replace('\n\n','\n').replace('\n','; '), delim)

