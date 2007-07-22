dnl  AMD K8 mpn_add_n/mpn_sub_n -- mpn add or subtract.

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


ifdef(`OPERATION_add_n', `
	define(M4_inst,        adcq)
	define(M4_function_n,  mpn_add_n)
	define(M4_function_nc, mpn_add_nc)
	define(M4_description, add)
',`ifdef(`OPERATION_sub_n', `
	define(M4_inst,        sbbq)
	define(M4_function_n,  mpn_sub_n)
	define(M4_function_nc, mpn_sub_nc)
	define(M4_description, subtract)
',`m4_error(`Need OPERATION_add_n or OPERATION_sub_n
')')')

MULFUNC_PROLOGUE(mpn_add_n mpn_add_nc mpn_sub_n mpn_sub_nc)


C mp_limb_t M4_function_n (mp_ptr dst, mp_srcptr src1, mp_srcptr src2,
C                         mp_size_t size);
C mp_limb_t M4_function_nc (mp_ptr dst, mp_srcptr src1, mp_srcptr src2,
C	                   mp_size_t size, mp_limb_t carry);
C
C Calculate src1,size M4_description src2,size, and store the result in
C dst,size.  The return value is the carry bit from the top of the result (1
C or 0).
C
C The _nc version accepts 1 or 0 for an initial carry into the low limb of
C the calculation.  Note values other than 1 or 0 here will lead to garbage
C results.
C
C This code runs at 1.64 cycles/limb, which is probably the best possible
C with plain integer operations.  Each limb is 2 loads and 1 store, and in
C one cycle the K7 can do two loads, or a load and a store, leading to 1.5
C c/l.

dnl  Must have UNROLL_THRESHOLD >= 2, since the unrolled loop can't handle 1.
ifdef(`PIC',`
deflit(UNROLL_THRESHOLD, 8)
',`
deflit(UNROLL_THRESHOLD, 8)
')

C  Input:	dst	%rdi
C		src1	%rsi   -> %rbx
C		src2	%rdx
C		size	%rcx
C		carry	%r8


define(PARAM_SIZE, %r8)
define(SAVE_RBP, %r9)
define(SAVE_RBX, %r10)

	TEXT
	ALIGN(32)

PROLOGUE(M4_function_nc)
	movq	%r8, %rax
	jmp	L(start)
EPILOGUE()

PROLOGUE(M4_function_n)

	xorq	%rax, %rax	C carry
L(start):
	cmpq	$UNROLL_THRESHOLD, %rcx

	jae	L(unroll)

	leaq	(%rsi,%rcx,8), %rsi
	leaq	(%rdx,%rcx,8), %rdx

	leaq	(%rdi,%rcx,8), %rdi
	negq	%rcx
	shrq	%rax

L(simple):
	C rax	scratch
	C rbx
	C rcx	counter
	C rdx	src2
	C rsi   src1
	C rdi	dst
	C rbp

	movq	(%rsi,%rcx,8), %rax
	M4_inst	(%rdx,%rcx,8), %rax
	movq	%rax, (%rdi,%rcx,8)
	incq	%rcx
	jnz	L(simple)

	movq	$0, %rax
	setc	%al

	ret

C -----------------------------------------------------------------------------
L(unroll):
	movq	%rbp, SAVE_RBP
	movq	%rcx, PARAM_SIZE   C for later use
	movq	%rbx, SAVE_RBX
	andq	$-2, %rcx		C size low bit masked out
	movq	%rsi, %rbx

	andq	$1, PARAM_SIZE		C size low bit kept

	movq 	%rdi, %rbp

	movq	%rcx, %rdi
	decq	%rcx

	shrq	$UNROLL_LOG2, %rcx
	negq	%rdi

	andq	$UNROLL_MASK, %rdi
	movq	%rdi, %rsi
	movq	%rdi, %r11
	shlq	$2, %r11
	shlq	$3, %rsi

ifdef(`PIC',`
	call	L(pic_calc)
L(here):
',`
	leaq	L(entry) (%rsi,%r11,1), %rsi	C 12 bytes per
')
	negq	%rdi
	shrq	%rax

	leaq	ifelse(UNROLL_BYTES,256,128) (%rbx,%rdi,8), %rbx
	leaq	ifelse(UNROLL_BYTES,256,128) (%rdx,%rdi,8), %rdx
	leaq	ifelse(UNROLL_BYTES,256,128) (%rbp,%rdi,8), %rdi

	jmp	*%rsi


ifdef(`PIC',`
L(pic_calc):
	C See mpn/x86/README about old gas bugs
	leaq	(%rsi,%r11,1), %rsi
	addq	$L(entry)-L(here), %rsi
	addq	(%rsp), %rsi
	ret
')


C -----------------------------------------------------------------------------
	ALIGN(32)
L(top):
	C rax	zero
	C rbx	src1
	C rcx	counter
	C rdx	src2
	C rsi	scratch (was computed jump)
	C rdi	dst
	C rbp	scratch

	leaq	UNROLL_BYTES(%rdx), %rdx

L(entry):
deflit(CHUNK_COUNT, 2)
forloop(i, 0, UNROLL_COUNT/CHUNK_COUNT-1, `
	deflit(`disp0', eval(i*CHUNK_COUNT*8 ifelse(UNROLL_BYTES,256,-128)))
	deflit(`disp1', eval(disp0 + 8))

Zdisp(	movq,	disp0,(%rbx), %rsi)
	movq	disp1(%rbx), %rbp
Zdisp(	M4_inst,disp0,(%rdx), %rsi)
Zdisp(	movq,	%rsi, disp0,(%rdi))
	M4_inst	disp1(%rdx), %rbp
	movq	%rbp, disp1(%rdi)
')

	decq	%rcx
	leaq	UNROLL_BYTES(%rbx), %rbx
	leaq	UNROLL_BYTES(%rdi), %rdi
	jns	L(top)


	movq	PARAM_SIZE, %rsi
	movq	$0, %rax

	decq	%rsi
	js	L(even)

	movq	(%rbx), %rcx
	M4_inst	UNROLL_BYTES(%rdx), %rcx
	movq	%rcx, (%rdi)
L(even):

	setc	%al

	movq	SAVE_RBP, %rbp
	movq	SAVE_RBX, %rbx

	ret


EPILOGUE()
