dnl  x86_64 mpn_sub_n -- Subtract two limb vectors of the same length > 0 and
dnl  store difference in a third limb vector.

dnl  Copyright 2004 Free Software Foundation, Inc.

dnl  This file is part of the GNU MP Library.

dnl  The GNU MP Library is free software; you can redistribute it and/or modify
dnl  it under the terms of the GNU Lesser General Public License as published
dnl  by the Free Software Foundation; either version 2.1 of the License, or (at
dnl  your option) any later version.

dnl  The GNU MP Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl  License for more details.

dnl  You should have received a copy of the GNU Lesser General Public License
dnl  along with the GNU MP Library; see the file COPYING.LIB.  If not, write
dnl  to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
dnl  Boston, MA 02110-1301, USA.

dnl **********************************************************************
dnl Actually, this file was created by Jason Martin <martin@math.jmu.edu>
dnl in an attempt to get GMP to run well on the Woodcrest CPU (aka Xeon).
dnl The GMP developers do not maintain this code and should not be
dnl bother with questions about it.  If you find errors in it, please
dnl let me know!
dnl **********************************************************************

dnl
dnl This is just a check to see if we are in my code testing sandbox
dnl or if we are actually in the GMP source tree
dnl
ifdef(`__JWM_Test_Code__',`
include(`./config.m4')
define(`MPN_PREFIX',`jwm_mpn_')',`
include(`../config.m4')')

dnl
dnl This is just a little test to see if the lahf and sahf
dnl instructions are available.  These instructions allow
dnl us to quickly save the EFLAGS register into AH and
dnl restore the EFLAGS from AH.  However, the earliest 64 bit
dnl CPUs didn't support this function, so the GNU assembler
dnl doesn't allow the lahf and sahf operands on 64 bit machines.
dnl To get around this, we check to see if the instructions are
dnl available.  If they are, then we use hand assembled bytes.
dnl If they aren't available then we fall back to using the
dnl setc and bt instructions which are slightly slower.
dnl
define(`LAHF_SAHF_AVAIL',esyscmd(`./lahf_sahf_test.sh'))
ifelse(LAHF_SAHF_AVAIL,`Yes',`
define(`save_CF_to_reg_a',`.byte	0x9f')
define(`get_CF_from_reg_a',`.byte	0x9e')',`
define(`save_CF_to_reg_a',`setc	%al')
define(`get_CF_from_reg_a',`bt	`$'0x0,%rax')')

C         cycles/limb
C Hammer:     2.5 (for 1024 limbs)
C Woodcrest:  2.6 (for 1024 limbs)

C INPUT PARAMETERS
C rp	rdi
C up	rsi
C vp	rdx
C n	rcx

ASM_START()
PROLOGUE(mpn_sub_n)
	pushq	%rbp		C Save off callee-save registers
	pushq	%rbx
	pushq	%r12
	pushq	%r13
	pushq	%r14
	pushq	%r15

	xor	%r15,%r15		C r15 will be our index, so
					C I'll call it i here after
	save_CF_to_reg_a				C Save CF

	mov	%rcx,%r9
	sub	$4,%r9			C r9 = n-(i+4)

	ALIGN(4)			C aligning for loop
L_mpn_sub_n_main_loop:
	C The goal of our main unrolled loop is to keep all the
	C execution units as busy as possible.  Since
	C there are three ALUs, we try to perform three
	C adds at a time.  Of course, we will have the
	C borrow dependency, so there is at least one
	C clock cycle between each sbb.  However, we'll
	C try to keep the other execution units busy
	C with loads and stores at the same time so that
	C our net throughput is close to one sbb per clock
	C cycle.  Hopefully this function will have asymptotic
	C behavior of taking 3*n clock cycles where n is the
	C number of limbs to sub.
	C
	C Note that I'm using FOUR sbbs at a time, this is just
	C because I wanted to use up all available registers since
	C I'm hoping the out-of-order and loop-pipeline logic in
	C the Xeon will help us out.

	C See if we are still looping
	jle	L_mpn_sub_n_loop_done

	get_CF_from_reg_a			C recover CF

	C Load inputs into rbx and r8
	C sub with borrow, and put result in r8
	C then store r8 to output.
	movq	(%rdx,%r15,8),%rbx
	movq	(%rsi,%r15,8),%r8
	sbbq	%rbx,%r8
	movq	%r8,(%rdi,%r15,8)

	C Load inputs into r9 and r10
	C sub with borrow, and put result in r10
	C then store r10 to output.
	movq	8(%rdx,%r15,8),%r9
	movq	8(%rsi,%r15,8),%r10
	sbbq	%r9,%r10
	movq	%r10,8(%rdi,%r15,8)

	C Load inputs into r11 and r12
	C sub with borrow, and put result in r12
	C then store r12 to output.
	movq	16(%rdx,%r15,8),%r11
	movq	16(%rsi,%r15,8),%r12
	sbbq	%r11,%r12
	movq	%r12,16(%rdi,%r15,8)

	C Load inputs into r13 and r14
	C sub with borrow, and put result in r14
	C then store r14 to output.
	movq	24(%rdx,%r15,8),%r13
	movq	24(%rsi,%r15,8),%r14
	sbbq	%r13,%r14
	movq	%r14,24(%rdi,%r15,8)

	save_CF_to_reg_a			C save CF

	mov	%r15,%r10
	add	$8,%r10
	add	$4,%r15		C increment by 4.

	mov	%rcx,%r9
	sub	%r10,%r9	C r9 = n-(i+4)
	jmp	L_mpn_sub_n_main_loop

L_mpn_sub_n_loop_done:
	mov	%rcx,%r15	C
	sub	%r9,%r15	C r15 = n-(n-(i+4))=i+4
	sub	$4,%r15		C r15 = i
	cmp	%rcx,%r15
L_mpn_sub_n_post_loop:
	je	L_mpn_sub_n_exit
	get_CF_from_reg_a			C recover CF

	C Load inputs into rbx and r8
	C sub with borrow, and put result in r8
	C then store r8 to output.
	movq	(%rdx,%r15,8),%rbx
	movq	(%rsi,%r15,8),%r8
	sbbq	%rbx,%r8
	movq	%r8,(%rdi,%r15,8)
	save_CF_to_reg_a			C save CF
	add	$1,%r15
	cmp	%rcx,%r15
	jmp	L_mpn_sub_n_post_loop


L_mpn_sub_n_exit:
	xor	%rcx,%rcx
	get_CF_from_reg_a			C recover the CF
	mov	%rcx,%rax	C Clears rax without affecting carry flag
	adc	%rax,%rax	C returns carry status.

	popq	%r15		C restore callee-save registers
	popq	%r14
	popq	%r13
	popq	%r12
	popq	%rbx
	popq	%rbp
	ret
EPILOGUE()
