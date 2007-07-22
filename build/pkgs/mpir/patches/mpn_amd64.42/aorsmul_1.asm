dnl  AMD K8 mpn_addmul_1/mpn_submul_1 -- add or subtract mpn multiple.

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

dnl 32 is the max.
deflit(UNROLL_COUNT, 16)

ifdef(`OPERATION_addmul_1',`
	define(M4_inst,        addq)
	define(M4_function_1,  mpn_addmul_1)
	define(M4_function_1c, mpn_addmul_1c)
	define(M4_description, add it to)
	define(M4_desc_retval, carry)
',`ifdef(`OPERATION_submul_1',`
	define(M4_inst,        subq)
	define(M4_function_1,  mpn_submul_1)
	define(M4_function_1c, mpn_submul_1c)
	define(M4_description, subtract it from)
	define(M4_desc_retval, borrow)
',`m4_error(`Need OPERATION_addmul_1 or OPERATION_submul_1
')')')

MULFUNC_PROLOGUE(mpn_addmul_1 mpn_addmul_1c mpn_submul_1 mpn_submul_1c)


C mp_limb_t M4_function_1 (mp_ptr dst, mp_srcptr src, mp_size_t size,
C                            mp_limb_t mult);
C mp_limb_t M4_function_1c (mp_ptr dst, mp_srcptr src, mp_size_t size,
C                             mp_limb_t mult, mp_limb_t carry);
C
C Calculate src,size multiplied by mult and M4_description dst,size.
C Return the M4_desc_retval limb from the top of the result.

ifdef(`PIC',`
deflit(UNROLL_THRESHOLD, 9)
',`
deflit(UNROLL_THRESHOLD, 6)
')

deflit(`FRAME',0)

define(PARAM_MULTIPLIER, %r9)
define(SAVE_RBX,        %r10)
define(SAVE_RBP,        %r11)


C  Input: dst         %rdi
C         src         %rsi
C         size        %rdx
C         multiplier  %rcx
C         carry       %r8

C  To be saved/preserved: %rbx, %rbp, %r12-15

	TEXT
	ALIGN(32)
PROLOGUE(M4_function_1)
	movq	%rcx, PARAM_MULTIPLIER
	movq	%rbx, SAVE_RBX
	movq	%rbp, SAVE_RBP

	movq	%rsi, %rax
	xorq	%rcx, %rcx

	decq	%rdx
	jnz	L(start_1)

	movq	(%rax), %rax
	movq	%rdi, %rcx

	mulq	PARAM_MULTIPLIER

	M4_inst	%rax, (%rcx)
	adcq	$0, %rdx
	movq	%rdx, %rax

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(M4_function_1c)
        movq    %rcx, PARAM_MULTIPLIER
        movq    %rbx, SAVE_RBX
        movq    %rbp, SAVE_RBP

	movq	%rsi, %rax

	decq	%rdx
	jnz	L(more_than_one_limb)

	movq	(%rax), %rax
	movq	%rdi, %rcx

	mulq	PARAM_MULTIPLIER

	addq	%r8, %rax

	adcq	$0, %rdx
	M4_inst	%rax, (%rcx)

	adcq	$0, %rdx
	movq	%rdx, %rax

	ret


	C offset ??? Should align ?
	ALIGN(16)
L(more_than_one_limb):
	movq	%r8, %rcx

	ALIGN(32)
L(start_1):
	C rax	src
	C rcx	initial carry
	C rdx	size-1
	C rsi   src
	C rdi   dst

	movq	%rdx, %rbx	C size-1

	cmpq	$UNROLL_THRESHOLD, %rdx

	movq	PARAM_MULTIPLIER, %rbp

	movq	(%rsi), %rax	C src low limb
	ja	L(unroll)


	C simple loop

	leaq	8(%rsi,%rbx,8), %rsi	C point one limb past last
	leaq	(%rdi,%rbx,8), %rdi	C point at last limb
	negq	%rbx

	C The movl to load the next source limb is done well ahead of the
	C mul.  This is necessary for full speed, and leads to one limb
	C handled separately at the end.

L(simple):
	C rax	src limb
	C rbx	loop counter
	C rcx	carry limb
	C rdx	scratch
	C rsi	src
	C rdi	dst
	C rbp	multiplier

	mulq	%rbp

	addq	%rax, %rcx
	adcq	$0, %rdx

	M4_inst	%rcx, (%rdi,%rbx,8)
	movq	(%rsi,%rbx,8), %rax
	adcq	$0, %rdx

	incq	%rbx
	movq	%rdx, %rcx
	jnz	L(simple)


	mulq	%rbp

	movq	SAVE_RBX, %rbx
	movq	SAVE_RBP, %rbp

	addq	%rax, %rcx
	adcq	$0, %rdx

	M4_inst	%rcx, (%rdi)
	adcq	$0, %rdx

	movq	%rdx, %rax
	ret

C -----------------------------------------------------------------------------
	ALIGN(32)
L(unroll):
	C rax	src low limb
	C rbx	size-1
	C rcx	carry
	C rdx	size-1
	C rsi	src
	C rdi	dst
	C rbp	multiplier

dnl  overlapping with parameters no longer needed
define(VAR_COUNTER,%r9)
define(VAR_JUMP,   %r8)

	subq	$2, %rbx	C (size-2)-1
	decq	%rdx		C size-2

	shrq	$UNROLL_LOG2, %rbx
	negq	%rdx

	movq	%rbx, VAR_COUNTER
	andq	$UNROLL_MASK, %rdx

	movq	%rdx, %r8
	movq	%rdx, %rbx
	shlq	$3, %r8
	shlq	$4, %rdx

ifdef(`PIC',`
	callq	L(pic_calc)
L(here):
',`
	leaq	L(entry) (%rdx, %r8, 1), %rdx
')
	negq	%rbx
	movq	%rdx, VAR_JUMP

	mulq	%rbp

	addq	%rax, %rcx	C initial carry, becomes low carry
	adcq	$0, %rdx
	testb	$1, %bl

	movq	8(%rsi), %rax	C src second limb
	leaq	ifelse(UNROLL_BYTES,256,128+) 16(%rsi,%rbx,8), %rsi
	leaq	ifelse(UNROLL_BYTES,256,128)   (%rdi,%rbx,8), %rdi

	movq	%rdx, %rbx	C high carry
	cmovnz	%rcx, %rbx	C high,low carry other way around
	cmovnz	%rdx, %rcx

	jmp	*VAR_JUMP


ifdef(`PIC',`
L(pic_calc):
	C See mpn/x86/README about old gas bugs
	leaq	(%rdx,%r8,1), %rdx
	addq	$L(entry)-L(here), %rdx
	addq	(%rsp), %rdx
	ret
')


C -----------------------------------------------------------------------------
C This code uses a "two carry limbs" scheme.  At the top of the loop the
C carries are ebx=lo, ecx=hi, then they swap for each limb processed.  For
C the computed jump an odd size means they start one way around, an even
C size the other.  Either way one limb is handled separately at the start of
C the loop.
C
C The positioning of the movl to load the next source limb is important.
C Moving it after the adcl with a view to avoiding a separate mul at the end
C of the loop slows the code down.
C [ This comment was for k7: for amd64, this is no longer the case! ]

	ALIGN(32)
L(top):
	C rax	src limb
	C rbx	carry high
	C rcx	carry low
	C rdx	scratch
	C rsi	src+8
	C rdi	dst
	C rbp	multiplier
	C
	C VAR_COUNTER  loop counter
	C
	C 24 bytes each limb

L(entry):
deflit(CHUNK_COUNT,2)
forloop(`i', 0, UNROLL_COUNT/CHUNK_COUNT-1, `
	deflit(`disp0', eval(i*CHUNK_COUNT*8 ifelse(UNROLL_BYTES,256,-128)))
	deflit(`disp1', eval(disp0 + 8))

	mulq	%rbp

Zdisp(	M4_inst,%rcx, disp0,(%rdi))
	movq	$0, %rcx

	adcq	%rax, %rbx

	adcq	%rdx, %rcx


Zdisp(	movq,	disp0,(%rsi), %rax)
	mulq	%rbp

	M4_inst	%rbx, disp1(%rdi)
	movq	$0, %rbx

	adcq	%rax, %rcx

	adcq	%rdx, %rbx
	movq	disp1(%rsi), %rax
')     C  ... END OF forloop

	decq	VAR_COUNTER
	leaq	UNROLL_BYTES(%rsi), %rsi
	leaq	UNROLL_BYTES(%rdi), %rdi

	jns	L(top)


	C rax	src limb
	C rbx	carry high
	C rcx	carry low
	C rdx
	C rsi
	C rdi	dst (points at second last limb)
	C rbp	multiplier
deflit(`disp0', ifelse(UNROLL_BYTES,256,-128))
deflit(`disp1', eval(disp0-0 + 8))

	mulq	%rbp

	M4_inst	%rcx, disp0(%rdi)
	movq	SAVE_RBP, %rbp

	adcq	%rbx, %rax
	movq	SAVE_RBX, %rbx

	adcq	$0, %rdx
	M4_inst	%rax, disp1(%rdi)

	adcq	$0, %rdx

	movq	%rdx, %rax
	ret


EPILOGUE()
