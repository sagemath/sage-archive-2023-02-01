from __future__ import with_statement

import ncadoctest
import sage.misc.randstate as randstate

OrigDocTestRunner = ncadoctest.DocTestRunner
class SageDocTestRunner(OrigDocTestRunner):
    def __init__(self, checker=None, verbose=None, optionflags=0):
        optionflags |= ncadoctest.NORMALIZE_WHITESPACE
        optionflags |= ncadoctest.ELLIPSIS
        OrigDocTestRunner.__init__(self, checker=checker, verbose=verbose, optionflags=optionflags)
        self._collect_timeit_stats = True
        self._timeit_stats = {}
        self._reset_random_seed = True
        self._random_seed = randstate.seed(0)

    def run(self, test, compileflags=None, out=None, clear_globs=True):
        r = OrigDocTestRunner.run(self, test, compileflags=compileflags, out=out, clear_globs=clear_globs)
        if self._collect_timeit_stats:
            pass # could save timeit stats here
        return r

    def run_one_test(self, test, compileflags, out):
        if self._reset_random_seed:
            randstate.set_random_seed(long(0))
        OrigDocTestRunner.run_one_test(self, test, compileflags, out)

    def run_one_example(self, test, example, filename, compileflags):
        if self._collect_timeit_stats:
            with self._random_seed:
                from sage.misc.sage_timeit import sage_timeit
                key = (example.source, example)
                try:
                    self._timeit_stats[key] = sage_timeit(example.source, test.globs)
                except Exception, e:
                    self._timeit_stats[key] = e
        # otherwise, just run the example
        OrigDocTestRunner.run_one_example(self, test, example, filename, compileflags)

    def save_timeit_stats_to_file_named(self, output_filename):
        if self._collect_timeit_stats:
            from sage.structure.sage_object import save
            save(self._timeit_stats, filename=output_filename)

ncadoctest.DocTestRunner = SageDocTestRunner

def testmod_returning_runner(m=None, name=None, globs=None, verbose=None,
                             report=True, optionflags=0, extraglobs=None,
                             raise_on_error=False, exclude_empty=False,
                             runner=None):
    return ncadoctest.testmod_returning_runner(m=m, name=name, globs=globs, verbose=verbose,
                                               report=report, optionflags=optionflags, extraglobs=extraglobs,
                                               raise_on_error=raise_on_error, exclude_empty=exclude_empty,
                                               runner=runner)
