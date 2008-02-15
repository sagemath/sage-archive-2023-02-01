#!/usr/bin/perl
# Script to answer queries of the form "is this available in Debian"
@use_debian = ('bzip2', 'atlas', 'blas', 'clisp', 'cvxopt', 'cython', 'ecm', 'f2c', 'fortran', 'freetype', 'gap', 'gd', 'gdmodule', 'gmp', 'gnutls', 'gsl', 'ipython', 'lapack', 'libgcrypt', 'libgpg_error', 'libpng', 'matplotlib', 'maxima', 'mercurial', 'moin', 'mpfi', 'mpfr', 'networkx', 'numpy', 'opencdk', 'pari', 'pexpect', 'pycrypto', 'python', 'python_gnutls', 'quaddouble', 'r', 'readline','scipy','scons','sqlite','sympy','termcap','twisted','weave','zlib','zodb3');
@our_debian = ('cddlib', 'eclib', 'flint', 'flintqs', 'libfplll', 'genus2reduction','gfan','givaro','iml','lcalc','libm4ri','linbox','ntl','palp','singular','symmetrica','sympow');
@ignore = (@use_debian,@our_debian);

while(<>) {
    $arg = $_;
    print "arg: $arg\n";
    foreach (@ignore) {
	if ($arg =~ m/^$_-/) {
	    $arg =~ s/.spkg$//;
	    print "$arg";
	    last;
	}
    }
}
