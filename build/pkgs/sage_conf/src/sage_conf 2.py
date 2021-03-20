# build/pkgs/sage_conf/src/sage_conf.py.  Generated from sage_conf.py.in by configure.

VERSION = "9.3.beta3"

MAXIMA = "/Users/charalampos/Desktop/Courses/PersonalProject/sage_mirror/local/bin/maxima"

ARB_LIBRARY = "arb"

SAGE_NAUTY_BINS_PREFIX = ""

# Colon-separated list of pkg-config modules to search for cblas functionality.
# We hard-code it here as cblas because configure (build/pkgs/openblas/spkg-configure.m4)
# always provides cblas.pc, if necessary by creating a facade pc file for a system BLAS.
CBLAS_PC_MODULES = "cblas"

# Used in sage.repl.ipython_kernel.install
JSMOL_DIR   = "/Users/charalampos/Desktop/Courses/PersonalProject/sage_mirror/local/share/jsmol"
MATHJAX_DIR = "/Users/charalampos/Desktop/Courses/PersonalProject/sage_mirror/local/share/mathjax"
THREEJS_DIR = "/Users/charalampos/Desktop/Courses/PersonalProject/sage_mirror/local/share/threejs"

# The following must not be used during build to determine source or installation
# location of sagelib.  See comments in SAGE_ROOT/src/Makefile.in
SAGE_LOCAL = "/Users/charalampos/Desktop/Courses/PersonalProject/sage_mirror/local"
SAGE_ROOT = "/Users/charalampos/Desktop/Courses/PersonalProject/sage_mirror"

# Entry point 'sage-config'.  It does not depend on any packages.

def _main():
    from argparse import ArgumentParser
    from sys import exit, stdout
    parser = ArgumentParser()
    parser.add_argument('--version', help="show version", action="version",
                       version='%(prog)s ' + VERSION)
    parser.add_argument("VARIABLE", nargs='?', help="output the value of VARIABLE")
    args = parser.parse_args()
    d = globals()
    if args.VARIABLE:
        stdout.write('{}\n'.format(d[args.VARIABLE]))
    else:
        for k, v in d.items():
            if not k.startswith('_'):
                stdout.write('{}={}\n'.format(k, v))

if __name__ == "__main__":
    _main()
