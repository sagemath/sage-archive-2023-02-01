dnl  AMD K8 mpn_mul_basecase -- multiply two mpn numbers.

dnl  This file is just an adaptation of similar file in the k7 directory.
dnl  Adapted by P. Gaudry in April 2005.
dnl  Here is the copyright of the original k7 version:

dnl  Copyright 1999, 2000, 2001, 2002 Free Software Foundation, Inc.
dnl
dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or
dnl  modify it under the terms of the GNU Lesser General Public License as
dnl  published by the Free Software Foundation; either version 2.1 of the
dnl  License, or (at your option) any later version.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful,
dnl  but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl  Lesser General Public License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public
dnl  License along with the GNU MP Library; see the file COPYING.LIB.  If
dnl  not, write to the Free Software Foundation, Inc., 59 Temple Place -
dnl  Suite 330, Boston, MA 02111-1307, USA.

include(`../config.m4')


deflit(UNROLL_COUNT, 16)


C void mpn_mul_basecase (mp_ptr wp,
C                        mp_srcptr xp, mp_size_t xsize,
C                        mp_srcptr yp, mp_size_t ysize);
C
C Calculate xp,xsize multiplied by yp,ysize, storing the result in
C wp,xsize+ysize.
C
C This routine is essentially the same as mpn/generic/mul_basecase.c, but
C it's faster because it does most of the mpn_addmul_1() startup
C calculations only once.  The saving is 15-25% on typical sizes coming from
C the Karatsuba multiply code.

ifdef(`PIC',`
deflit(UNROLL_THRESHOLD, 5)
',`
deflit(UNROLL_THRESHOLD, 5)
')

define(param_ysize,%r8)
define(param_yp,   %rcx)
define(param_xsize,%rdx)
define(param_xp,   %rsi)
define(param_wp,   %rdi)


	TEXT
	ALIGN(32)
PROLOGUE(mpn_mul_basecase)
deflit(`FRAME',0)

	cmpq	$2, param_xsize
	ja	L(xsize_more_than_two)
	je	L(two_by_something)

	C one limb by one limb

	movq	(param_yp), %rax	C rax <- y0
	mulq	(param_xp)              C [ax:dx] <- y0*x0

	movq	%rax, (param_wp)
	movq	%rdx, 8(param_wp)
	ret


C -----------------------------------------------------------------------------
L(two_by_something):
C xsize = 2, hence rdx is free for usage
deflit(`FRAME',0)
	decq	param_ysize		C YSIZE--

	movq	(param_yp), %r9		C r9 <- y0
	movq	(param_xp), %rax	C rax <- x0
	jnz	L(two_by_two)


	C two limbs by one limb

	mulq	%r9			C  [ax:dx] <- x0*y0

	movq	%rax, (param_wp)    	C  w0 <- low_prod
	movq	8(param_xp), %rax   	C  rax <- x1   (rsi is now free)
	movq	%rdx, %rsi		C  rsi <- carry

	mulq	%r9			C  [ax:dx] <- x1*y0

	addq	%rax, %rsi		C  rsi <- ax+carry   ( --> carry_flag)
	movq	%rsi, 8(param_wp)   	C  w1 <- rsi

	adcq	$0, %rdx		C  dx <- dx+carry
	movq	%rdx, 16(param_wp)  	C  w2 <- dx

	ret


C -----------------------------------------------------------------------------
	ALIGN(16)
L(two_by_two):
	C rax	x0			r8
	C rbx	**untouched**		r9      y0
	C rcx	yp                      r10-11
	C rdx
	C rsi	xp
	C rdi   wp
	C rbp

	mulq	%r9			C [ax:dx]  <- x0*y0

	movq	%rdx, %r10		C r10 <- carry for w1

	movq	%rax, (param_wp)	C w0 <- ax
	movq	8(param_xp), %rax	C ax <- x1

	mulq	%r9			C [ax:dx]  <- x1*y0

	addq	%rax, %r10		C r10 <- r10 + ax  for w1

	adcq	$0, %rdx		C dx <- carry for w2
	movq	8(param_yp), %rcx	C cx <- y1
	movq	%r10, 8(param_wp)	C w1 <- r10

	movq	8(param_xp), %rax	C ax <- x1
	movq	%rdx, %r10		C carry, for w2

	mulq	%rcx			C [ax:dx] <- x1*y1

	addq	%rax, %r10      	C r10 <- for w2

	adcq	$0, %rdx		C for w3
	movq	(param_xp), %rax	C x0

	movq	%rdx, %rsi		C carry, for w3

	mulq	%rcx			C x0*y1

	addq	%rax, 8(param_wp)	C w1 += ax
	adcq	%rdx, %r10		C for w2
	movq	%r10, 16(param_wp)	C w2 <- r10

	adcq	$0, %rsi
	movq	%rsi, 24(param_wp)	C w3 <- carry in rsi

	ret


C -----------------------------------------------------------------------------
	ALIGN(16)
L(xsize_more_than_two):

C The first limb of yp is processed with a simple mpn_mul_1 style loop
C inline.  Unrolling this doesn't seem worthwhile since it's only run once
C (whereas the addmul below is run ysize-1 many times).  A call to the
C actual mpn_mul_1 will be slowed down by the call and parameter pushing and
C popping, and doesn't seem likely to be worthwhile on the typical 13-26
C limb operations the Karatsuba code calls here with.

	C rax			r8	ysize
	C rbx
	C rcx	yp
	C rdx	xsize
	C rsi	xp
	C rdi	wp
	C rbp


define(PARAM_YSIZE,%r8)    C already there
define(PARAM_YP,   %r9)    C  init : %rcx
define(PARAM_XSIZE,%r10)   C  init : %rdx
define(PARAM_XP,   %r11)   C  init : %rsi
define(PARAM_WP,   %r12)   C  init : %rdi     r12 should be saved!


dnl  FRAME doesn't carry on from previous, no pushes yet here
defframe(`SAVE_RBX',-8)
defframe(`SAVE_R12',-16)
defframe(`SAVE_RBP',-24)
deflit(`FRAME',0)

	subq	$24, %rsp
deflit(`FRAME',24)

	movq	%rbx, SAVE_RBX
	movq	%r12, SAVE_R12
	movq	%rbp, SAVE_RBP

	movq	%rcx, PARAM_YP
	movq	%rdx, PARAM_XSIZE
	movq	%rsi, PARAM_XP
	movq	%rdi, PARAM_WP

	movq	(PARAM_YP), %rbp
	movq	PARAM_XSIZE, %rcx

	xorq	%rbx, %rbx
	leaq	(PARAM_XP,PARAM_XSIZE,8), %rsi	C xp end

	leaq	(PARAM_WP,PARAM_XSIZE,8), %rdi	C wp end of mul1
	negq	%rcx


L(mul1):
	C rax	scratch
	C rbx	carry
	C rcx	counter, negative
	C rdx	scratch
	C rsi	xp end
	C rdi	wp end of mul1
	C rbp	multiplier

	movq	(%rsi,%rcx,8), %rax

	mulq	%rbp

	addq	%rbx, %rax
	movq	%rax, (%rdi,%rcx,8)
	movq	$0, %rbx

	adcq	%rdx, %rbx
	incq	%rcx
	jnz	L(mul1)


	movq	PARAM_YSIZE, %rdx
	movq	PARAM_XSIZE, %rcx

	movq	%rbx, (%rdi)		C final carry
	decq	%rdx

	jnz	L(ysize_more_than_one)


	movq	SAVE_RBX, %rbx
	movq	SAVE_RBP, %rbp
	movq	SAVE_R12, %r12
	addq	$FRAME, %rsp

	ret


L(ysize_more_than_one):
	cmpq	$UNROLL_THRESHOLD, %rcx
	movq	PARAM_YP, %rax

	jae	L(unroll)


C -----------------------------------------------------------------------------
	C simple addmul looping
	C
	C rax	yp
	C rbx
	C rcx	xsize
	C rdx	ysize-1
	C rsi	xp end
	C rdi	wp end of mul1
	C rbp

	leaq	8(%rax,%rdx,8), %rbp	C yp end
	negq	%rcx
	negq	%rdx

	movq	(%rsi,%rcx,8), %rax	C xp low limb
	movq	%rdx, PARAM_YSIZE	C -(ysize-1)
	incq	%rcx

	xorq	%rbx, %rbx		C initial carry
	movq	%rcx, PARAM_XSIZE	C -(xsize-1)
	movq	%rbp, PARAM_YP

	movq	(%rbp,%rdx,8), %rbp	C yp second lowest limb - multiplier
	jmp	L(simple_outer_entry)


	C this is offset ????  Align ?
L(simple_outer_top):
	C rbp	ysize counter, negative

	movq	PARAM_YP, %rdx
	movq	PARAM_XSIZE, %rcx	C -(xsize-1)
	xorq	%rbx, %rbx		C carry

	movq	%rbp, PARAM_YSIZE
	addq	$8, %rdi		C next position in wp

	movq	(%rdx,%rbp,8), %rbp	C yp limb - multiplier
	movq	-8(%rsi,%rcx,8), %rax	C xp low limb


L(simple_outer_entry):

L(simple_inner):
	C rax	xp limb
	C rbx	carry limb
	C rcx	loop counter (negative)
	C rdx	scratch
	C rsi	xp end
	C rdi	wp end
	C rbp	multiplier

	mulq	%rbp

	addq	%rax, %rbx
	adcq	$0, %rdx

	addq	%rbx, (%rdi,%rcx,8)
	movq	(%rsi,%rcx,8), %rax
	adcq	$0, %rdx

	incq	%rcx
	movq	%rdx, %rbx
	jnz	L(simple_inner)


	mulq	%rbp

	movq	PARAM_YSIZE, %rbp
	addq	%rax, %rbx

	adcq	$0, %rdx
	addq	%rbx, (%rdi)

	adcq	$0, %rdx
	incq	%rbp

	movq	%rdx, 8(%rdi)
	jnz	L(simple_outer_top)


	movq	SAVE_RBX, %rbx
	movq	SAVE_RBP, %rbp
	movq	SAVE_R12, %r12
	addq	$FRAME, %rsp

	ret



C -----------------------------------------------------------------------------
C
C The unrolled loop is the same as in mpn_addmul_1(), see that code for some
C comments.
C
C VAR_ADJUST is the negative of how many limbs the leals in the inner loop
C increment xp and wp.  This is used to adjust back xp and wp, and rshifted
C to given an initial VAR_COUNTER at the top of the outer loop.
C
C VAR_COUNTER is for the unrolled loop, running from VAR_ADJUST/UNROLL_COUNT
C up to -1, inclusive.
C
C VAR_JMP is the computed jump into the unrolled loop.
C
C VAR_XP_LOW is the least significant limb of xp, which is needed at the
C start of the unrolled loop.
C
C PARAM_YSIZE is the outer loop counter, going from -(ysize-1) up to -1,
C inclusive.
C
C PARAM_YP is offset appropriately so that the PARAM_YSIZE counter can be
C added to give the location of the next limb of yp, which is the multiplier
C in the unrolled loop.
C
C The trick with VAR_ADJUST means it's only necessary to do one fetch in the
C outer loop to take care of xp, wp and the inner loop counter.

defframe(VAR_COUNTER,  -32)
defframe(VAR_ADJUST,   -40)
defframe(VAR_XP_LOW,   -48)
deflit(VAR_EXTRA_SPACE, 24)


L(unroll):
	C rax	yp
	C rbx
	C rcx	xsize
	C rdx	ysize-1
	C rsi	xp end
	C rdi	wp end of mul1
	C rbp

	movq	PARAM_XP, %rsi          C from here, PARAM_XP not used
	movq	8(%rax), %rbp		C multiplier (yp second limb)
	leaq	8(%rax,%rdx,8), %rax	C yp adjust for ysize indexing

	movq	PARAM_WP, %rdi
	movq	%rax, PARAM_YP
	negq	%rdx

		C  From here, only PARAM_YP and PARAM_YSIZE are used
		C  Hence r10, r11, r12 are free for use

	movq	%rdx, PARAM_YSIZE
	leaq	UNROLL_COUNT-2(%rcx), %rbx	C (xsize-1)+UNROLL_COUNT-1
	decq	%rcx				C xsize-1

	movq	(%rsi), %rax		C xp low limb
	andq	$-UNROLL_MASK-1, %rbx
	negq	%rcx

	subq	$VAR_EXTRA_SPACE, %rsp
deflit(`FRAME',24+VAR_EXTRA_SPACE)
	negq	%rbx
	andq	$UNROLL_MASK, %rcx
	movq	%rcx, %r12		C for later parity test

	movq	%rbx, VAR_ADJUST
	movq	%rcx, %rdx
	movq	%rcx, %r10
	shlq	$4, %rcx
	shlq	$3, %r10

	sarq	$UNROLL_LOG2, %rbx

	C 24=16+8 code bytes per limb
ifdef(`PIC',`
	callq	L(pic_calc)
L(unroll_here):
',`
	leaq	L(unroll_entry) (%rcx,%r10,1), %rcx
')
	negq	%rdx

	movq	%rax, VAR_XP_LOW
	movq	%rcx, PARAM_XP		C PARAM_XP used for VAR_JUMP
	leaq	8(%rdi,%rdx,8), %rdi	C wp and xp, adjust for unrolling,
	leaq	8(%rsi,%rdx,8), %rsi	C  and start at second limb
	jmp	L(unroll_outer_entry)


ifdef(`PIC',`
L(pic_calc):
	C See mpn/x86/README about old gas bugs
	leaq	(%rcx,%r10,1), %rcx
	addq	$L(unroll_entry)-L(unroll_here), %rcx
	addq	(%rsp), %rcx
	ret
')


C --------------------------------------------------------------------------
	ALIGN(32)
L(unroll_outer_top):
	C ebp	ysize counter, negative

	movq	VAR_ADJUST, %rbx
	movq	PARAM_YP, %rdx

	movq	VAR_XP_LOW, %rax
	movq	%rbp, PARAM_YSIZE	C store incremented ysize counter

	leaq	8(%rdi,%rbx,8), %rdi
	leaq	(%rsi,%rbx,8), %rsi
	sarq	$UNROLL_LOG2, %rbx

	movq	(%rdx,%rbp,8), %rbp	C yp next multiplier
	movq	PARAM_XP, %rcx

L(unroll_outer_entry):
	mulq	%rbp

	movq	%r12, %rcx
	testb	$1, %cl		C and clear carry bit
	movq	%rbx, VAR_COUNTER
	movq	$0, %rbx

	movq	$0, %rcx
	cmovz	%rax, %rcx	C eax into low carry, zero into high carry limb
	cmovnz	%rax, %rbx

	C Extra fetch of VAR_JMP is bad, but registers are tight
	C TODO: we have more registers, now!!!!
	jmp	*PARAM_XP


C -----------------------------------------------------------------------------
	ALIGN(32)
L(unroll_top):
	C rax	xp limb
	C rbx	carry high
	C rcx	carry low
	C rdx	scratch
	C rsi	xp+8
	C rdi	wp
	C rbp	yp multiplier limb
	C
	C VAR_COUNTER  loop counter, negative
	C
	C 24 bytes each limb

L(unroll_entry):

deflit(CHUNK_COUNT,2)
forloop(`i', 0, UNROLL_COUNT/CHUNK_COUNT-1, `
	deflit(`disp0', eval(i*CHUNK_COUNT*8 ifelse(UNROLL_BYTES,256,-128)))
	deflit(`disp1', eval(disp0 + 8))

	adcq	%rdx, %rbx
Zdisp(	movq,	disp0,(%rsi), %rax)

	mulq	%rbp

Zdisp(	addq,	%rcx, disp0,(%rdi))
	movq	$0, %rcx

	adcq	%rax, %rbx


	adcq	%rdx, %rcx
	movq	disp1(%rsi), %rax

	mulq	%rbp

	addq	%rbx, disp1(%rdi)
	movq	$0, %rbx

	adcq	%rax, %rcx
')


	incq	VAR_COUNTER
	leaq	UNROLL_BYTES(%rsi), %rsi
	leaq	UNROLL_BYTES(%rdi), %rdi

	jnz	L(unroll_top)


	C rax
	C rbx	zero
	C rcx	low
	C rdx	high
	C rsi
	C rdi	wp, pointing at second last limb)
	C rbp
	C
	C carry flag to be added to high

deflit(`disp0', ifelse(UNROLL_BYTES,256,-128))
deflit(`disp1', eval(disp0-0 + 8))

	movq	PARAM_YSIZE, %rbp
	adcq	$0, %rdx
	addq	%rcx, disp0(%rdi)

	adcq	$0, %rdx
	incq	%rbp

	movq	%rdx, disp1(%rdi)
	jnz	L(unroll_outer_top)


	movq	SAVE_RBP, %rbp
	movq	SAVE_R12, %r12
	movq	SAVE_RBX, %rbx
	addq	$FRAME, %rsp

	ret

EPILOGUE()
