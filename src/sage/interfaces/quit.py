"""nodoctest"""

import os

expect_objects = []

def expect_quitall(verbose=False):
    for P in expect_objects:
        R = P()
        if not R is None:
            try:
                R.quit(verbose=verbose)
            except RuntimeError, msg:
                if verbose:
                    print msg




