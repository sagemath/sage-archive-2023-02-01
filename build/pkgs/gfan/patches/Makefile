ALL: default

# GPROFFLAG = -pg
# GPROFFLAG = -g
#PREFIX?=/usr/local
#BINDIR=$(PREFIX)/bin

RELEASEDIR = gfan0.4plus
MAIN       = gfan

GCATSPATH   = ./


ifeq ($(gmppath),)
GMP_LINKOPTIONS = -lgmp
GMP_INCLUDEOPTIONS =
else
GMP_LINKOPTIONS = $(gmppath)/lib/libgmp.a
GMP_INCLUDEOPTIONS = -I $(gmppath)/include
endif

ifeq ($(cddpath),)
CDD_LINKOPTIONS = -lcddgmp
CDD_INCLUDEOPTIONS =
else
CDD_LINKOPTIONS = $(cddpath)/lib/libcddgmp.a
CDD_INCLUDEOPTIONS = -I $(cddpath)/include
endif

ifeq ($(soplex),)
SOPLEX_PATH =
SOPLEX_LINKOPTIONS =
SOPLEX_INCLUDEOPTIONS =
SOPLEX_OBJECTS =
else
SOPLEX_PATH = $(HOME)/math/software/soplex-1.3.2
# SOPLEX_LINKOPTIONS = -lz $(SOPLEX_PATH)/lib/libsoplex.linux.x86.gnu.opt.a
SOPLEX_LINKOPTIONS = -lz $(SOPLEX_PATH)/lib/libsoplex.darwin.x86.gnu.opt.a
SOPLEX_INCLUDEOPTIONS = -I $(SOPLEX_PATH)/src
SOPLEX_OBJECTS = lp_soplexcdd.o
endif

ifeq ($(singular),)
SINGULAR_PATH =
SINGULAR_LINKOPTIONS =
SINGULAR_INCLUDEOPTIONS =
SINGULAR_OBJECTS =
else
SINGULAR_PATH = $(HOME)/math/software/Singular-3-1-0
SINGULAR_LINKOPTIONS =  -L$(SINGULAR_PATH)/Singular -lsingular -lncurses -lreadline
SINGULAR_INCLUDEOPTIONS = -I $(SINGULAR_PATH)/kernel -I $(SINGULAR_PATH)/omalloc
SINGULAR_OBJECTS = singular.o
endif



ADDITIONALLINKOPTIONS = $(CDD_LINKOPTIONS) $(GMP_LINKOPTIONS) $(SOPLEX_LINKOPTIONS) $(SINGULAR_LINKOPTIONS)
ADDITIONALINCLUDEOPTIONS = $(CDD_INCLUDEOPTIONS) $(GMP_INCLUDEOPTIONS) $(SOPLEX_INCLUDEOPTIONS) $(SINGULAR_INCLUDEOPTIONS)


MKDIR=mkdir -p


SHELL       = /bin/sh
ARCH        = LINUX
CC          = gcc
CLINKER     = $(CC)
CCC         = g++
CCLINKER    = $(CCC)
#OPTFLAGS    = -O2 -DGMPRATIONAL -DNDEBUG
OPTFLAGS    = -O2 -DGMPRATIONAL

CFLAGS	  = -I$(SAGE_LOCAL)/include $(OPTFLAGS) $(GPROFFLAG) $(ADDITIONALINCLUDEOPTIONS) #-pedantic
CCFLAGS	  = $(CFLAGS)
FFLAGS	  = $(OPTFLAGS)

CATSOBJECTS =	$(GCATSPATH)lp_cdd.o \
		$(SOPLEX_OBJECTS) \
		$(SINGULAR_OBJECTS) \
		$(GCATSPATH)parser.o \
		$(GCATSPATH)field.o \
		$(GCATSPATH)monomial.o \
		$(GCATSPATH)printer.o	\
		$(GCATSPATH)polynomial.o \
		$(GCATSPATH)termorder.o \
		$(GCATSPATH)term.o \
		$(GCATSPATH)vektor.o \
		$(GCATSPATH)division.o \
		$(GCATSPATH)buchberger.o \
		$(GCATSPATH)wallideal.o \
		$(GCATSPATH)lp.o \
		$(GCATSPATH)enumeration.o \
		$(GCATSPATH)ep_standard.o \
		$(GCATSPATH)ep_xfig.o \
		$(GCATSPATH)reversesearch.o \
		$(GCATSPATH)application.o \
		$(GCATSPATH)timer.o \
		$(GCATSPATH)renderer.o \
		$(GCATSPATH)field_rationals.o \
		$(GCATSPATH)symmetry.o \
		$(GCATSPATH)breadthfirstsearch.o \
		$(GCATSPATH)genericwalk.o \
		$(GCATSPATH)minkowskisum.o \
		$(GCATSPATH)newtonpolytope.o \
		$(GCATSPATH)tropical.o \
		$(GCATSPATH)dimension.o \
		$(GCATSPATH)bergman.o \
		$(GCATSPATH)subspace.o \
		$(GCATSPATH)polyhedralcone.o \
		$(GCATSPATH)gfanapplication.o \
		$(GCATSPATH)polyhedralfan.o \
		$(GCATSPATH)tropical2.o \
		$(GCATSPATH)gmpallocator.o \
		$(GCATSPATH)field_zmodpz.o \
		binomial.o \
		matrix.o \
		latticeideal.o \
		scarf.o \
		xfig.o \
		halfopencone.o \
		lll.o \
		multiplicity.o \
		substitute.o \
		polymakefile.o \
		saturation.o \
		determinant.o \
		polynomialring.o \
		log.o \
		tropicalbasis.o \
		symmetriccomplex.o \
		linalg.o \
		minors.o \
		continuedfractions.o \
		triangulation.o \
		minkowskidual.o \
		regularsubdivision.o \
		fieldlp.o \
		field_rationalfunctions.o \
		tropical_weildivisor.o \
		intsinpolytope.o\
		lattice.o \
		graph.o \
		restrictedautoreduction.o \
		tropicaltraverse.o \
		groebnerengine.o \
		ge_gfan.o \
		nbody.o \
		codimoneconnectedness.o \
	#	restrictedgfan.o \

APPOBJECTS = app_main.o \
		app_buchberger.o \
		app_construction.o \
		app_doesidealcontain.o \
		app_edwinsconjecture.o \
		app_facets.o \
		app_fvector.o \
		app_groebnercone.o \
		app_groupfacetbinomials.o \
		app_homogeneityspace.o \
		app_homogenize.o \
		app_initialforms.o \
		app_interactive.o \
		app_isgroebnerbasis.o \
		app_ismarkedgroebnerbasis.o \
		app_krulldimension.o \
		app_leadingterms.o \
		app_liststandardmonomials.o \
		app_markpolynomialset.o \
		app_multiplymatrix.o \
		app_polynomialsetunion.o \
		app_render.o \
		app_rendernewtonpolytope.o \
		app_renderstaircase.o \
		app_stats.o \
		app_substitute.o \
		app_supportindices.o \
		app_tolatex.o \
		app_transposematrix.o \
		app_tropical.o \
		app_tropicalbasis.o \
		app_tropicalintersection.o \
		app_tropicalstartingcone.o \
		app_tropicaltraverse.o \
		app_walk.o \
		app_weightvector.o \
		app_xfigconstruction.o \
		app_scarfisgeneric.o \
		app_scarfvisualize.o \
		app_scarfcomplex.o \
		app_sturmsequence.o \
		app_latticeideal.o \
		app_latticetest.o \
		app_lll.o \
		app_tropicalmultiplicity.o \
		app_idealintersection.o \
		app_test.o \
		app_saturation.o \
		app_idealproduct.o \
		app_representatives.o \
		app_tropicallifting.o \
		app_add.o \
		app_topolyhedralfan.o \
		app_tropicalbruteforce.o \
		app_secondaryfan.o \
		app_composepermutations.o \
		app_minors.o \
		app_tropicalrank.o \
		app_minkowski.o \
		app_polytopetopolynomial.o \
		app_triangulate.o \
		app_tropicallinearspace.o \
		app_grassmanndata.o \
		app_combinerays.o \
		app_regularsubdivision.o \
		app_berndssuggestion.o \
		app_grassmanndata2.o \
		app_grassmanndata3.o \
		app_lpsolve.o \
		app_tropicalweildivisor.o \
		app_lattice.o \
		app_intsinpolytope.o\
		app_tropicalevaluation.o \
		app_smalessixth.o \
		app_smalessixth2.o \
		app_moeckel.o \
		app_nbody.o \
		app_spolynomial.o \
		app_link.o \
		app_normalfancleanup.o \
		app_istriangulation.o \
		app_tropicalfunction.o \
		app_volume.o \
		app_isconnected.o \
		app_checkridges.o \
		app_tropicalhypersurface.o \
		app_product.o \
		app_commonrefinement.o \
#		app_isrefinement.o \    # needs to be fixed so that it compiles with gcc version 2.96 (legolas.imf.au.dk)


EXECS	  = $(MAIN)

# Define suffixes to make the program compile on legolas.imf.au.dk :
.SUFFIXES: .o .cpp .c

OBJECTS = $(CATSOBJECTS) $(APPOBJECTS)

all: $(MAIN)

$(BINDIR): $(PREFIX)
	$(MKDIR) $(BINDIR)

$(PREFIX):
	$(MKDIR) $(PREFIX)

default: $(OBJECTS) $(ADDITIONALOBJECTS) $(EXECS)

$(MAIN): $(OBJECTS)
	$(CCLINKER) $(OBJECTS) $(ADDITIONALLINKOPTIONS) $(GPROFFLAG) -o $(MAIN)

release:
	rm -f -r $(RELEASEDIR)/*
	mkdir -p $(RELEASEDIR)
	cp *.cpp $(RELEASEDIR)
	cp *.h $(RELEASEDIR)
	cp Makefile $(RELEASEDIR)
	cp README $(RELEASEDIR)
	cp LICENSE $(RELEASEDIR)
	cp COPYING $(RELEASEDIR)
	cp Doxyfile $(RELEASEDIR)
	mkdir -p $(RELEASEDIR)/examples/
#General examples:
	cp examples/2x2of2x3 $(RELEASEDIR)/examples/
	cp examples/2x2of2x4 $(RELEASEDIR)/examples/
	cp examples/2x2of3x3 $(RELEASEDIR)/examples/
	cp examples/2x2of4x4 $(RELEASEDIR)/examples/
	cp examples/4x4of4x5 $(RELEASEDIR)/examples/
	cp examples/4x4of5x5 $(RELEASEDIR)/examples/
	cp examples/6x6-subPfaffians $(RELEASEDIR)/examples/
	cp examples/cyclic4 $(RELEASEDIR)/examples/
	cp examples/linhyper5_2 $(RELEASEDIR)/examples/
	cp examples/linhyper5_2.cone $(RELEASEDIR)/examples/
	cp examples/pablo $(RELEASEDIR)/examples/
	cp examples/symmetryTest $(RELEASEDIR)/examples/
#Examples in Groebner fan paper:
	cp examples/examplePaper $(RELEASEDIR)/examples/
	cp examples/sturmfels3.9 $(RELEASEDIR)/examples/
	cp examples/3x3of3x4 $(RELEASEDIR)/examples/
	cp examples/3x3of3x5 $(RELEASEDIR)/examples/
	cp examples/3x3of4x4 $(RELEASEDIR)/examples/
#	cp examples/3x3of4x4sym $(RELEASEDIR)/examples/
	cp examples/grassmann2_5 $(RELEASEDIR)/examples/
	cp examples/cyclic5 $(RELEASEDIR)/examples/
#	cp examples/J4 $(RELEASEDIR)/examples/
#Examples useful for tropical computations:
#	cp examples/grassmann2_5 $(RELEASEDIR)/examples/
	cp examples/grassmann2_5.cone $(RELEASEDIR)/examples/
	cp examples/grassmann2_6 $(RELEASEDIR)/examples/
	cp examples/grassmann2_6.cone $(RELEASEDIR)/examples/
	cp examples/grassmann3_6 $(RELEASEDIR)/examples/
	cp examples/grassmann3_6.cone $(RELEASEDIR)/examples/
#Examples in tropical paper:
	cp examples/hankel3x3of4x4 $(RELEASEDIR)/examples/
	cp examples/hankel3x3of4x4.cone $(RELEASEDIR)/examples/
	cp examples/hankel3x3of4x5 $(RELEASEDIR)/examples/
	cp examples/hankel3x3of4x5.cone $(RELEASEDIR)/examples/
#	cp examples/hankel3x3of5x5 $(RELEASEDIR)/examples/
#	cp examples/hankel3x3of5x5.cone $(RELEASEDIR)/examples/
	cp examples/3x3of3x5.cone $(RELEASEDIR)/examples/
#	cp examples/3x3of4x4sym $(RELEASEDIR)/examples/
	cp examples/3x3of4x4sym.cone $(RELEASEDIR)/examples/
#	cp examples/3x3of5x5sym $(RELEASEDIR)/examples/
#	cp examples/3x3of5x5sym.cone $(RELEASEDIR)/examples/
	cp examples/commat2x2 $(RELEASEDIR)/examples/
	cp examples/commat2x2.cone $(RELEASEDIR)/examples/
#	cp examples/commat3x3 $(RELEASEDIR)/examples/
#	cp examples/commat3x3.cone $(RELEASEDIR)/examples/

	mkdir -p $(RELEASEDIR)/doc/
	cp doc/Makefile $(RELEASEDIR)/doc/
	cp doc/*.bib $(RELEASEDIR)/doc/
	cp doc/*.bbl $(RELEASEDIR)/doc/
	cp doc/*.tex $(RELEASEDIR)/doc/
	cp doc/*.dvi $(RELEASEDIR)/doc/
	cp doc/*.eps $(RELEASEDIR)/doc/
	cp doc/*.bst $(RELEASEDIR)/doc/
	cp doc/Makefile $(RELEASEDIR)/doc/
	mkdir -p $(RELEASEDIR)/homepage/
	cp homepage/*.png $(RELEASEDIR)/homepage/
	cp homepage/*.html $(RELEASEDIR)/homepage/
	cp homepage/Makefile $(RELEASEDIR)/homepage/
	mkdir -p $(RELEASEDIR)/homepage/presentation
	cp homepage/presentation/*.fig $(RELEASEDIR)/homepage/presentation/
	cp homepage/presentation/*.tex $(RELEASEDIR)/homepage/presentation/
	cp homepage/presentation/*.eps $(RELEASEDIR)/homepage/presentation/
	cp homepage/presentation/*.bib $(RELEASEDIR)/homepage/presentation/
	cp homepage/presentation/*.bbl $(RELEASEDIR)/homepage/presentation/
	cp homepage/presentation/*.ps $(RELEASEDIR)/homepage/presentation/
	cp homepage/presentation/*.pdf $(RELEASEDIR)/homepage/presentation/
	cp homepage/presentation/Makefile $(RELEASEDIR)/homepage/presentation/

	tar -c $(RELEASEDIR) > $(RELEASEDIR).tar
	gzip $(RELEASEDIR).tar

clean:
	/bin/rm -f *.o $(EXECS) $(MAIN)
install: $(BINDIR)
	cp $(EXECS) $(BINDIR)
#	cp $(EXECS) /usr/local/bin
#	./gfan installlinks --path $(BINDIR)/
	cd $(BINDIR) && ./gfan installlinks
.c.o:
	$(CC) $(CFLAGS) -c $<
.cc.o:
	$(CCC) -c $<
.cpp.o:
	$(CCC) $(CFLAGS) -c $<
.C.o:
	$(CCC) -c $<
# wget http://ftp.sunet.se/pub/gnu/gmp/gmp-4.2.2.tar.gz
# tar -xzvf gmp-4.2.2.tar.gz
# cd gmp-4.2.2
# ./configure --prefix=$HOME/gmp
# make
# make install
# make check
# cd..

# wget ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlib-094d.tar.gz  # THIS LINE DOES NOT WORK!
# tar -xzvf cddlib-094d.tar.gz
# cd cddlib-094d
# ./configure --prefix="$HOME/cddlib" CFLAGS="-I$HOME/gmp/include -L$HOME/gmp/lib"
