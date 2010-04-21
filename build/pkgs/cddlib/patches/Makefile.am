##############################################
# Site customizations:
#---------------------------------------------
# where is gmp?
#---------------------------------------------
# export gmpdir := $(PWD)/../gmp-4.1
export gmpdir := $(SAGE_LOCAL)
#---------------------------------------------
# where is gmp lib?
export gmplibdir = $(gmpdir)/lib
#---------------------------------------------
# where is gmp include?
#---------------------------------------------
# export gmpincludedir = $(gmpdir)
export gmpincludedir = $(gmpdir)/include
#---------------------------------------------
# End of site customization.
##############################################
ACLOCAL_AMFLAGS = -I m4
SUBDIRS      = lib-src src lib-src-gmp src-gmp
EXTRA_DIST   = examples/*.* examples-ine/*.* examples-ine3d/*.* examples-ml/*.* examples-ext/*.*
