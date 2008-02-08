ALL: default

# GPROFFLAG = -pg

RELEASEDIR = gfan0.2.2
MAIN       = gfan

GCATSPATH   = ./

ifeq ($(cddpath),)
ADDITIONALLINKOPTIONS = -lcddgmp -lgmp
#ADDITIONALLINKOPTIONS = -lgmp -lcddgmp
ADDITIONALINCLUDEOPTIONS =
else
#ADDITIONALLINKOPTIONS = -lgmp $(cddpath)lib-src-gmp/libcddgmp.a
ADDITIONALLINKOPTIONS = $(cddpath)lib-src-gmp/libcddgmp.a -lgmp
ADDITIONALINCLUDEOPTIONS = -I $(cddpath)lib-src -I $(cddpath)lib-src-gmp
endif

SHELL       = /bin/sh
ARCH        = LINUX
CC          = gcc
CLINKER     = $(CC)
CCC         = g++
CCLINKER    = $(CCC)
OPTFLAGS    = -O2 -DGMPRATIONAL

CFLAGS	  = -I$(SAGE_LOCAL)/include -L$(SAGE_LOCAL)/lib $(OPTFLAGS) $(GPROFFLAG) $(ADDITIONALINCLUDEOPTIONS) #-pedantic
#CFLAGS	  = $(OPTFLAGS) $(GPROFFLAG) $(ADDITIONALINCLUDEOPTIONS) #-pedantic
CCFLAGS	  = $(CFLAGS)
FFLAGS	  = $(OPTFLAGS)

CATSOBJECTS =	$(GCATSPATH)lp_cdd.o \
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
		$(GCATSPATH)field_zmodpz.o \
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
		binomial.o \
		matrix.o \
		latticeideal.o \
		scarf.o \
		xfig.o \


# I removed this from the list below
#    app_construction.o \
# See SAGE.txt for more details.

APPOBJECTS = app_main.o \
		app_buchberger.o \
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
		app_isrefinement.o \
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
		app_latticeideal.o \

EXECS	  = $(MAIN)

# Define suffixes to make the program compile on legolas.imf.au.dk :
.SUFFIXES: .o .cpp .c

OBJECTS = $(CATSOBJECTS) $(APPOBJECTS)

all: $(MAIN)

default: $(OBJECTS) $(ADDITIONALOBJECTS) $(EXECS)

$(MAIN): $(OBJECTS)
	$(CCLINKER) $(CFLAGS) $(OBJECTS) $(ADDITIONALLINKOPTIONS) $(GPROFFLAG) -o $(MAIN)

release:
	rm $(RELEASEDIR)/* -f -r
	mkdir $(RELEASEDIR) -p
	cp *.cpp $(RELEASEDIR)
	cp *.h $(RELEASEDIR)
	cp Makefile $(RELEASEDIR)
	cp README $(RELEASEDIR)
	cp LICENSE $(RELEASEDIR)
	cp COPYING $(RELEASEDIR)
	mkdir $(RELEASEDIR)/examples/ -p
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

	mkdir $(RELEASEDIR)/doc/ -p
	cp doc/Makefile $(RELEASEDIR)/doc/
	cp doc/*.bib $(RELEASEDIR)/doc/
	cp doc/*.bbl $(RELEASEDIR)/doc/
	cp doc/*.tex $(RELEASEDIR)/doc/
	cp doc/*.dvi $(RELEASEDIR)/doc/
	cp doc/*.eps $(RELEASEDIR)/doc/
	cp doc/*.bst $(RELEASEDIR)/doc/
	cp doc/Makefile $(RELEASEDIR)/doc/
	mkdir $(RELEASEDIR)/homepage/ -p
	cp homepage/*.png $(RELEASEDIR)/homepage/
	cp homepage/*.html $(RELEASEDIR)/homepage/
	cp homepage/Makefile $(RELEASEDIR)/homepage/
	mkdir $(RELEASEDIR)/homepage/presentation -p
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
install:
	cp $(EXECS) /usr/local/bin
	./gfan installlinks /usr/local/bin/
.c.o:
	$(CC) $(CFLAGS) -c $<
.cc.o:
	$(CCC) -c $<
.cpp.o:
	$(CCC) $(CFLAGS) -c $<
.C.o:
	$(CCC) -c $<
