from lazy_attribute import lazy_attribute, lazy_class_attribute
from lazy_import import lazy_import

from misc import (alarm, cancel_alarm,
                  ellipsis_range, ellipsis_iter, srange, xsrange, sxrange, getitem,
                  cputime, verbose, set_verbose, set_verbose_files,
                  get_verbose_files, unset_verbose_files, get_verbose,
                  version, banner, add, union, uniq, powerset, subsets,
                  exists, forall, is_iterator,
                  random_sublist, mul, walltime, generic_cmp,
                  repr_lincomb,
                  pad_zeros, attrcall,
                  SAGE_DB, SAGE_TMP,
                  is_32_bit, is_64_bit, newton_method_sizes, compose,
                  self_compose, nest)

from temporary_file import tmp_dir, tmp_filename

from misc_c import prod, running_total, balanced_sum

from dev_tools import runsnake, import_statements

from html import html

from table import table

from sage_timeit_class import timeit

from edit_module import edit, set_edit_template

from flatten import flatten

from map_threaded import map_threaded

from session import load_session, save_session, show_identifiers

from remote_file import get_remote_file

from profiler import Profiler

from mrange import xmrange, mrange, xmrange_iter, mrange_iter, cartesian_product_iterator

from fpickle import pickle_function, unpickle_function

# deprecated
#from bug import bug

from dist import install_scripts

from package import install_package, is_package_installed, standard_packages, optional_packages, experimental_packages, upgrade

from pager import pager

from sagedoc import (search_src, search_def, search_doc, browse_sage_doc,
                     tutorial, reference, manual, developer, constructions,
                     python_help, help)

from classgraph import class_graph

from reset import reset, restore

from getusage import top, get_memory_usage

from log import log_html, log_dvi, log_text

from mathml import mathml

from defaults import set_default_variable_name

from preparser import preparse, implicit_multiplication, BackslashOperator

lazy_import('sage.misc.attached_files', [
        'attach', 'detach', 'attached_files', 'load_attach_path',
        'reset_load_attach_path', 'load_attach_mode'])

from interpreter import preparser

from sage_eval import sage_eval, sageobj

from sage_input import sage_input

from cython import cython_lambda, cython_create_local_so
from cython_c import cython
pyrex = cython # synonym -- for now
sagex = cython # synonym -- for now

from prun import prun

from persist import save, load, dumps, loads, db, db_save

from func_persist import func_persist

from functional import (additive_order,
                        sqrt as numerical_sqrt,
                        base_ring,
                        base_field,
                        basis,
                        category,
                        charpoly,
                        characteristic_polynomial,
                        coerce,
                        cyclotomic_polynomial,
                        decomposition,
                        denominator,
                        det,
                        dimension,
                        dim,
                        discriminant,
                        disc,
                        eta,
                        exp,
                        factor,
                        fcp,
                        gen,
                        gens,
                        hecke_operator,
                        ideal,
                        image,
                        integral, integrate,
                        integral_closure,
                        interval,
                        xinterval,
                        is_commutative,
                        is_even,
                        is_integrally_closed,
                        is_field,
                        is_odd,
                        kernel,
                        krull_dimension,
                        lift,
                        minimal_polynomial,
                        minpoly,
                        multiplicative_order,
                        ngens,
                        norm,
                        numerator,
                        numerical_approx,
                        n, N,
                        objgens,
                        objgen,
                        one,
                        order,
                        rank,
                        regulator,
                        round,
                        quotient,
                        quo,
                        show,
                        isqrt,
                        squarefree_part,
                        symbolic_sum as sum,
                        transpose,
                        zero,
                        log as log_b,
                        parent)


from latex import LatexExpr, latex, view, pretty_print, pretty_print_default

from trace import trace

from constant_function import ConstantFunction

from cachefunc import CachedFunction, cached_function, cached_method, cached_in_parent_method, disk_cached_function

from abstract_method import abstract_method

from randstate import seed, set_random_seed, initial_seed, current_randstate

from prandom import *

from sage_unittest import TestSuite

from explain_pickle import explain_pickle, unpickle_newobj, unpickle_global, unpickle_build, unpickle_instantiate, unpickle_persistent, unpickle_extension, unpickle_appends

from bitset import FrozenBitset, Bitset

from decorators import specialize, sage_wraps, infix_operator

from unknown import Unknown

##########################################################################
def benchmark(n=-1):
    """
    Run a well-chosen range of Sage commands and record the time it
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


import messaging

from ascii_art import ascii_art
