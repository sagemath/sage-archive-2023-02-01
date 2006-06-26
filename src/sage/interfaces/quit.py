"""nodoctest"""

import os

expect_objects = []

def expect_quitall(verbose=False):
    for P in expect_objects:
        R = P()
        if not R is None:
            if verbose:
                print "Quitting %s"%R
            try:
                R.quit()
            except RuntimeError, msg:
                if verbose:
                    print msg




