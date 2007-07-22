#
# This is a Makefile to generate test code for the optimizations
# I'm attempting to make to GMP for Core2 Macs.
#
#

M4=m4
CC=gcc
CFLAGS="-m64"
OBJS= mpn_test.o performance_mon.o orig_add_n.o add_n.o\
sub_n.o orig_sub_n.o addmul_1.o orig_addmul_1.o\
orig_submul_1.o submul_1.o orig_mul_1.o mul_1.o

all:	mpn_test

%.o : %.s
	$(CC) $(CFLAGS) -c $< -o $@

%.s : %.asm
	$(M4) -D__JWM_Test_Code__ $< > $@

mpn_test: $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean :
	rm -f mpn_test
	rm -f *.o
	rm -f *~
