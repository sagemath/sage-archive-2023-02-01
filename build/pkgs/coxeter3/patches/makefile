# sources contains a list of the source files (i.e., the .cpp files)
sources := $(patsubst %.cpp,%.cpp,$(wildcard *.cpp))
# there is one .o file for each .cpp file
objects := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
# this is used to automatically generate the dependencies
dependencies := $(patsubst %.cpp,%.d,$(wildcard *.cpp))

globals = globals.h

pflags = -c $(includedirs) -pg -O
oflags = -c $(includedirs) -O -Wall
gflags = -c $(includedirs) -g

cflags = $(gflags) # the default setting
cflags = -c $(includedirs) -O2 -fPIC -g -DFAST -DALLTRUE

ifdef optimize
	NDEBUG = true
	cflags = $(oflags)
endif

ifdef profile
	cflags = $(pflags)
endif

cc = g++

all: coxeter executable

coxeter: $(objects)
	$(cc) -shared -Wl,-soname,libcoxeter3.so -o libcoxeter3.so $(objects) -lc

executable: $(objects)
	$(cc) -o coxeter $(objects)

clean:
	rm -f $(objects)

%.o:%.cpp
	$(cc) $(cflags) $*.cpp

# dependencies --- these were generated automatically by make depend on my
# system; they are explicitly copied for portability. Only local dependencies
# are considered. If you add new #include directives you should add the
# corresponding dependencies here; the best way if your compiler supports the
# -MM option is probably to simply say make depend > tmp, then paste in the
# contents of tmp in lieu of the dependencies listed here.

%.d:%.cpp
	@$(cc) -MM $*.cpp
depend: $(dependencies)

affine.o: affine.cpp affine.h globals.h coxgroup.h coxtypes.h io.h list.h \
  memory.h constants.h list.hpp error.h files.h hecke.h interface.h \
  automata.h bits.h minroots.h dotval.h graph.h type.h transducer.h \
  schubert.h stack.h stack.hpp hecke.hpp polynomials.h vector.h \
  vector.hpp polynomials.hpp invkl.h klsupport.h search.h search.hpp kl.h \
  uneqkl.h wgraph.h files.hpp cells.h
automata.o: automata.cpp automata.h globals.h bits.h list.h memory.h \
  constants.h list.hpp error.h io.h
bits.o: bits.cpp bits.h globals.h list.h memory.h constants.h list.hpp \
  error.h io.h
cells.o: cells.cpp cells.h globals.h bits.h list.h memory.h constants.h \
  list.hpp error.h io.h kl.h coxtypes.h klsupport.h polynomials.h \
  vector.h vector.hpp polynomials.hpp schubert.h interface.h automata.h \
  minroots.h dotval.h graph.h type.h transducer.h stack.h stack.hpp \
  hecke.h hecke.hpp search.h search.hpp uneqkl.h wgraph.h
commands.o: commands.cpp commands.h globals.h coxgroup.h coxtypes.h io.h \
  list.h memory.h constants.h list.hpp error.h files.h hecke.h \
  interface.h automata.h bits.h minroots.h dotval.h graph.h type.h \
  transducer.h schubert.h stack.h stack.hpp hecke.hpp polynomials.h \
  vector.h vector.hpp polynomials.hpp invkl.h klsupport.h search.h \
  search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h dictionary.h \
  dictionary.hpp directories.h fcoxgroup.h help.h interactive.h special.h \
  typeA.h
constants.o: constants.cpp constants.h globals.h
coxgroup.o: coxgroup.cpp coxgroup.h globals.h coxtypes.h io.h list.h \
  memory.h constants.h list.hpp error.h files.h hecke.h interface.h \
  automata.h bits.h minroots.h dotval.h graph.h type.h transducer.h \
  schubert.h stack.h stack.hpp hecke.hpp polynomials.h vector.h \
  vector.hpp polynomials.hpp invkl.h klsupport.h search.h search.hpp kl.h \
  uneqkl.h wgraph.h files.hpp cells.h
coxtypes.o: coxtypes.cpp coxtypes.h globals.h io.h list.h memory.h \
  constants.h list.hpp error.h
error.o: error.cpp error.h globals.h coxtypes.h io.h list.h memory.h \
  constants.h list.hpp directories.h graph.h bits.h type.h interactive.h \
  interface.h automata.h minroots.h dotval.h transducer.h kl.h \
  klsupport.h polynomials.h vector.h vector.hpp polynomials.hpp \
  schubert.h stack.h stack.hpp hecke.h hecke.hpp search.h search.hpp \
  version.h
fcoxgroup.o: fcoxgroup.cpp fcoxgroup.h globals.h coxgroup.h coxtypes.h \
  io.h list.h memory.h constants.h list.hpp error.h files.h hecke.h \
  interface.h automata.h bits.h minroots.h dotval.h graph.h type.h \
  transducer.h schubert.h stack.h stack.hpp hecke.hpp polynomials.h \
  vector.h vector.hpp polynomials.hpp invkl.h klsupport.h search.h \
  search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h
files.o: files.cpp files.h globals.h hecke.h list.h memory.h constants.h \
  list.hpp error.h interface.h automata.h bits.h io.h coxtypes.h \
  minroots.h dotval.h graph.h type.h transducer.h schubert.h stack.h \
  stack.hpp hecke.hpp polynomials.h vector.h vector.hpp polynomials.hpp \
  invkl.h klsupport.h search.h search.hpp kl.h uneqkl.h wgraph.h \
  files.hpp cells.h directories.h posets.h version.h
general.o: general.cpp general.h globals.h coxgroup.h coxtypes.h io.h \
  list.h memory.h constants.h list.hpp error.h files.h hecke.h \
  interface.h automata.h bits.h minroots.h dotval.h graph.h type.h \
  transducer.h schubert.h stack.h stack.hpp hecke.hpp polynomials.h \
  vector.h vector.hpp polynomials.hpp invkl.h klsupport.h search.h \
  search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h
graph.o: graph.cpp graph.h globals.h list.h memory.h constants.h list.hpp \
  error.h bits.h io.h coxtypes.h type.h directories.h interactive.h \
  interface.h automata.h minroots.h dotval.h transducer.h
hecke.o: hecke.cpp
help.o: help.cpp help.h globals.h commands.h coxgroup.h coxtypes.h io.h \
  list.h memory.h constants.h list.hpp error.h files.h hecke.h \
  interface.h automata.h bits.h minroots.h dotval.h graph.h type.h \
  transducer.h schubert.h stack.h stack.hpp hecke.hpp polynomials.h \
  vector.h vector.hpp polynomials.hpp invkl.h klsupport.h search.h \
  search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h dictionary.h \
  dictionary.hpp directories.h
interactive.o: interactive.cpp interactive.h globals.h bits.h list.h \
  memory.h constants.h list.hpp error.h io.h coxtypes.h graph.h type.h \
  interface.h automata.h minroots.h dotval.h transducer.h affine.h \
  coxgroup.h files.h hecke.h schubert.h stack.h stack.hpp hecke.hpp \
  polynomials.h vector.h vector.hpp polynomials.hpp invkl.h klsupport.h \
  search.h search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h \
  directories.h fcoxgroup.h general.h typeA.h
interface.o: interface.cpp interface.h globals.h automata.h bits.h list.h \
  memory.h constants.h list.hpp error.h io.h coxtypes.h minroots.h \
  dotval.h graph.h type.h transducer.h
invkl.o: invkl.cpp invkl.h globals.h coxtypes.h io.h list.h memory.h \
  constants.h list.hpp error.h klsupport.h polynomials.h vector.h \
  vector.hpp polynomials.hpp schubert.h bits.h interface.h automata.h \
  minroots.h dotval.h graph.h type.h transducer.h stack.h stack.hpp \
  hecke.h hecke.hpp search.h search.hpp
io.o: io.cpp io.h globals.h list.h memory.h constants.h list.hpp error.h
kl.o: kl.cpp kl.h globals.h coxtypes.h io.h list.h memory.h constants.h \
  list.hpp error.h klsupport.h polynomials.h vector.h vector.hpp \
  polynomials.hpp schubert.h bits.h interface.h automata.h minroots.h \
  dotval.h graph.h type.h transducer.h stack.h stack.hpp hecke.h \
  hecke.hpp search.h search.hpp iterator.h
klsupport.o: klsupport.cpp klsupport.h globals.h coxtypes.h io.h list.h \
  memory.h constants.h list.hpp error.h polynomials.h vector.h vector.hpp \
  polynomials.hpp schubert.h bits.h interface.h automata.h minroots.h \
  dotval.h graph.h type.h transducer.h stack.h stack.hpp
main.o: main.cpp constants.h globals.h commands.h coxgroup.h coxtypes.h \
  io.h list.h memory.h list.hpp error.h files.h hecke.h interface.h \
  automata.h bits.h minroots.h dotval.h graph.h type.h transducer.h \
  schubert.h stack.h stack.hpp hecke.hpp polynomials.h vector.h \
  vector.hpp polynomials.hpp invkl.h klsupport.h search.h search.hpp kl.h \
  uneqkl.h wgraph.h files.hpp cells.h dictionary.h dictionary.hpp \
  version.h
memory.o: memory.cpp memory.h globals.h constants.h error.h
minroots.o: minroots.cpp minroots.h globals.h bits.h list.h memory.h \
  constants.h list.hpp error.h io.h coxtypes.h dotval.h graph.h type.h
polynomials.o: polynomials.cpp
posets.o: posets.cpp posets.h globals.h bits.h list.h memory.h \
  constants.h list.hpp error.h io.h wgraph.h interface.h automata.h \
  coxtypes.h minroots.h dotval.h graph.h type.h transducer.h
sage.o: sage.cpp sage.h globals.h coxtypes.h coxgroup.h list.h schubert.h
schubert.o: schubert.cpp schubert.h globals.h coxtypes.h io.h list.h \
  memory.h constants.h list.hpp error.h bits.h interface.h automata.h \
  minroots.h dotval.h graph.h type.h transducer.h stack.h stack.hpp
search.o: search.cpp
special.o: special.cpp special.h globals.h commands.h coxgroup.h \
  coxtypes.h io.h list.h memory.h constants.h list.hpp error.h files.h \
  hecke.h interface.h automata.h bits.h minroots.h dotval.h graph.h \
  type.h transducer.h schubert.h stack.h stack.hpp hecke.hpp \
  polynomials.h vector.h vector.hpp polynomials.hpp invkl.h klsupport.h \
  search.h search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h \
  dictionary.h dictionary.hpp directories.h interactive.h
stack.o: stack.cpp
transducer.o: transducer.cpp transducer.h globals.h coxtypes.h io.h \
  list.h memory.h constants.h list.hpp error.h graph.h bits.h type.h
type.o: type.cpp type.h globals.h io.h list.h memory.h constants.h \
  list.hpp error.h
typeA.o: typeA.cpp typeA.h globals.h fcoxgroup.h coxgroup.h coxtypes.h \
  io.h list.h memory.h constants.h list.hpp error.h files.h hecke.h \
  interface.h automata.h bits.h minroots.h dotval.h graph.h type.h \
  transducer.h schubert.h stack.h stack.hpp hecke.hpp polynomials.h \
  vector.h vector.hpp polynomials.hpp invkl.h klsupport.h search.h \
  search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h
uneqkl.o: uneqkl.cpp uneqkl.h globals.h coxtypes.h io.h list.h memory.h \
  constants.h list.hpp error.h hecke.h interface.h automata.h bits.h \
  minroots.h dotval.h graph.h type.h transducer.h schubert.h stack.h \
  stack.hpp hecke.hpp polynomials.h vector.h vector.hpp polynomials.hpp \
  klsupport.h search.h search.hpp interactive.h
vector.o: vector.cpp
wgraph.o: wgraph.cpp wgraph.h globals.h list.h memory.h constants.h \
  list.hpp error.h bits.h io.h interface.h automata.h coxtypes.h \
  minroots.h dotval.h graph.h type.h transducer.h stack.h stack.hpp
