lcalc: cmdline.c
	$(CXX) $(CFLAGS) -DINCLUDE_PARI   $(DEFINES)\
      -I$(SAGE_LOCAL)/include/pari -I$(SAGE_LOCAL)/include\
      -I ../include/ -L$(SAGE_LOCAL)/lib \
      cmdline.c \
      Lcommandline.cc Lcommandline_elliptic.cc Lcommandline_globals.cc \
      Lcommandline_misc.cc Lcommandline_numbertheory.cc \
      Lcommandline_twist.cc Lcommandline_values_zeros.cc \
      Lgamma.cc Lglobals.cc Lmisc.cc Lriemannsiegel.cc \
            -o lcalc -lpari -lmpfr -lgmpxx -lgmp
