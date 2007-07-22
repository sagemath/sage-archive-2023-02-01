dnl *********************************************************************
dnl  Intel64 mpn_submul_1 -- Multiply a limb vector with a limb and
dnl  subtract the result from a second limb vector.
dnl
dnl  Copyright 2006  Jason Worth Martin <jason.worth.martin@gmail.com>
dnl
dnl  This code is distributed under the terms of the Gnu Lesser
dnl  General Public License.
dnl
dnl  This code is part of a patch for the GNU MP Library (GMP).
dnl  However, this code is not maintained by the GMP group, so please
dnl  do not bother them with questions/bugs/etc. related to this code.
dnl  If you find a bug in this code, please contact Jason Worth Martin
dnl  at jason.worth.martin@gmail.com
dnl
dnl  This patch is free software; you can redistribute it and/or modify
dnl  it under the terms of the GNU Lesser General Public License as published
dnl  by the Free Software Foundation; either version 2.1 of the License, or (at
dnl  your option) any later version.
dnl
dnl  This patch is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl  License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public License
dnl  along with the GNU MP Library; see the file LICENSE.TXT.  If not, write
dnl  to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
dnl  Boston, MA 02110-1301, USA.
dnl
dnl *********************************************************************
dnl
dnl
dnl CREDITS
dnl
dnl This code is based largely on Pierrick Gaudry's excellent assembly
dnl support for the AMD64 architecture.  (Note that Intel64 and AMD64,
dnl while using the same instruction set, have very different
dnl microarchitectures.  So, this code performs very poorly on AMD64
dnl machines even though it is near-optimal on Intel64.)
dnl
dnl Roger Golliver works for Intel and provided insightful improvements
dnl particularly in using the "lea" instruction to perform additions
dnl and register-to-register moves.
dnl
dnl Eric Bainville has a brilliant exposition of optimizing arithmetic for
dnl AMD64 (http://www.bealto.it).  I adapted many of the ideas he
dnl describes to Intel64.
dnl
dnl Agner Fog is a demigod in the x86 world.  If you are reading assembly
dnl code files and you haven't heard of Agner Fog, then take a minute to
dnl look over his software optimization manuals (http://www.agner.org/).
dnl They are superb.
dnl


dnl With a 4-way unroll (UNROLL_EXPONENT = 2) the code has
dnl
dnl         	cycles/limb
dnl Hammer:           4.6
dnl Woodcrest:        4.6
dnl
dnl With increased unrolling, it appears to converge to 4 cycles/limb
dnl on Intel Core 2 machines.  I believe that this is optimal, however
dnl it requires such absurd unrolling that it becomes unusable for all
dnl but the largest inputs.  A 4-way unroll seems like a good balance
dnl to me because then commonly used input sizes (e.g. 1024bit Public
dnl keys) still benifit from the speed up.

dnl
dnl This is just a check to see if we are in my code testing sandbox
dnl or if we are actually in the GMP source tree
dnl
ifdef(`__JWM_Test_Code__',`
include(`./config.m4')
define(`MPN_PREFIX',`jwm_mpn_')',`
include(`../config.m4')')


dnl *********************************************************************
dnl *********************************************************************
dnl
dnl Here are the important m4 parameters for the code
dnl
dnl     BpL is Bytes per Limb (8 since this is 64bit code)
dnl
dnl	UNROLL_EXPONENT is the power of 2 for which we unroll the code.
dnl                     possible values are 1,2,...,8.  A reasonable
dnl                     value is probably 2.  If really large inputs
dnl                     are expected, then 4 is probably good.  Larger
dnl                     values are really only useful for flashy
dnl                     benchmarks and testing asymptotic behavior.
dnl
dnl     THRESHOLD is the minimum number of limbs needed before we bother
dnl               using the complicated loop.  A reasonable value is
dnl               2*UNROLL_SIZE + 6
dnl
dnl *********************************************************************
dnl *********************************************************************
define(`BpL',`8')
define(`UNROLL_EXPONENT',`2')
define(`UNROLL_SIZE',eval(2**UNROLL_EXPONENT))
define(`UNROLL_MASK',eval(UNROLL_SIZE-1))
define(`THRESHOLD',eval(2*UNROLL_SIZE+6))

dnl Here is a convenient Macro for addressing
dnl memory.  Entries of the form
dnl
dnl      ADDR(ptr,index,displacement)
dnl
dnl get converted to
dnl
dnl      displacement*BpL(ptr,index,BpL)
dnl
define(`ADDR',`eval(`$3'*BpL)($1,$2,BpL)')



C Register	Usage
C --------	-----
C rax		low word from mul
C rbx*
C rcx		s2limb
C rdx		high word from mul
C rsi		s1p
C rdi		rp
C rbp*		Base Pointer
C rsp*		Stack Pointer
C r8		A_x
C r9		A_y
C r10		A_z
C r11		B_x
C r12*		B_y
C r13*		B_z
C r14*		temp
C r15*		index
C
C * indicates that the register must be
C preserved for the caller.
define(`s2limb',`%rcx')
define(`s1p',`%rsi')
define(`rp',`%rdi')
define(`A_x',`%r8')
define(`A_y',`%r9')
define(`A_z',`%r10')
define(`B_x',`%r11')
define(`B_y',`%r12')
define(`B_z',`%r13')
define(`temp',`%r14')
define(`index',`%r15')


C INPUT PARAMETERS
C rp		rdi
C s1p		rsi
C n		rdx
C s2limb	rcx

ASM_START()
PROLOGUE(mpn_submul_1)
					C Compare the limb count
					C with the threshold value.
					C If the limb count is small
					C we just use the small loop,
					C otherwise we jump to the
					C more complicated loop.
	cmp	$THRESHOLD,%rdx
	jge	L_mpn_submul_1_main_loop_prep
	movq	%rdx, %r11
	leaq	(%rsi,%rdx,8), %rsi
	leaq	(%rdi,%rdx,8), %rdi
	negq	%r11
	xorq	%r8, %r8
	xorq	%r10, %r10
	jmp	L_mpn_submul_1_small_loop
	ALIGN(16)
L_mpn_submul_1_small_loop:
	movq	(%rsi,%r11,8), %rax
	mulq	%rcx
	addq	%r8, %rax
	adcq	%r10, %rdx
	subq	%rax, (%rdi,%r11,8)
	movq	%r10, %r8
	adcq	%rdx, %r8
	incq	%r11
	jne	L_mpn_submul_1_small_loop

	movq	%r8, %rax
	ret

L_mpn_submul_1_main_loop_prep:
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
				C If n is even, we need to do three
				C pre-multiplies, if n is odd we only
				C need to do two.
	movq	%rdx, temp
	movq	$0, index
	movq	$0, A_x
	movq	$0, A_y
	andq	$1, %rdx
	jnz	L_mpn_submul_1_odd_n

					C Case n is even
	movq	ADDR(s1p,index,0), %rax
	mulq	s2limb
	subq	%rax, ADDR(rp,index,0)
	adcq	%rdx, A_x
	addq	$1, index
					C At this point
					C  temp = n (even)
					C index = 1

L_mpn_submul_1_odd_n:
					C Now
					C temp = n
					C index = 1 if n even
					C       = 0 if n odd
					C
	movq	ADDR(s1p,index,0), %rax
	mulq	s2limb
	movq	ADDR(rp,index,0), A_z
	addq	%rax, A_x
	adcq	%rdx, A_y

	movq	ADDR(s1p,index,1), %rax
	mulq	s2limb
	movq	ADDR(rp,index,1), B_z
	movq	%rax, B_x
	movq	%rdx, B_y
	movq	ADDR(s1p,index,2), %rax

	addq	$3, index
	leaq	ADDR(s1p,temp,-1), s1p
	leaq	ADDR(rp,temp,-1), rp
	negq	temp
	addq	temp, index
				C At this point:
				C s1p   = address of last s1limb
				C rp    = address of last rplimb
				C temp  = -n
				C index = 4 - n if n even
				C       = 3 - n if n odd
				C
				C So, index is a (negative) even
				C number.
				C
				C *****************************************
				C ATTENTION:
				C
				C From here on, I will use array
				C indexing notation in the comments
				C because it is convenient.  So, I
				C will pretend that index is positive
				C because then a comment like
				C      B_z = rp[index-1]
				C is easier to read.
				C However, keep in mind that index is
				C actually a negative number indexing
				C back from the end of the array.
				C This is a common trick to remove one
				C compare operation from the main loop.
				C *****************************************

				C
				C Now we enter a spin-up loop the
				C will make sure that the index is
				C a multiple of UNROLL_SIZE before
				C going to our main unrolled loop.
	movq	index, temp
	negq	temp
	andq	$UNROLL_MASK, temp
	jz	L_mpn_submul_1_main_loop
	shrq	$1, temp
L_mpn_submul_1_main_loop_spin_up:
				C At this point we should have:
				C
				C A_x = low_mul[index-2] + carry_in
				C A_y = high_mul[index-2] + CF
				C A_z = rp[index-2]
				C
				C B_x = low_mul[index-1]
				C B_y = high_mul[index-1]
				C B_z = rp[index-1]
				C
				C rax = s1p[index]
	mulq	s2limb
	subq	A_x, A_z
	leaq	(%rax), A_x
	movq	ADDR(s1p,index,1), %rax
	adcq	A_y, B_x
	movq	A_z, ADDR(rp,index,-2)
	movq	ADDR(rp,index,0), A_z
	adcq	$0, B_y
	leaq	(%rdx), A_y
				C At this point we should have:
				C
				C B_x = low_mul[index-1] + carry_in
				C B_y = high_mul[index-1] + CF
				C B_z = rp[index-1]
				C
				C A_x = low_mul[index]
				C A_y = high_mul[index]
				C A_z = rp[index]
				C
				C rax = s1p[index+1]
	mulq	s2limb
	subq	B_x, B_z
	leaq	(%rax), B_x
	movq	ADDR(s1p,index,2), %rax
	adcq	B_y, A_x
	movq	B_z, ADDR(rp,index,-1)
	movq	ADDR(rp,index,1), B_z
	adcq	$0, A_y
	leaq	(%rdx), B_y

	addq	$2, index
	subq	$1, temp
	jnz	L_mpn_submul_1_main_loop_spin_up

	jmp	L_mpn_submul_1_main_loop
	ALIGN(16)
L_mpn_submul_1_main_loop:
				C The code here is really the same
				C logic as the spin-up loop.  It's
				C just been unrolled.
forloop(`unroll_index', 0, eval(UNROLL_SIZE/2 - 1),`
	mulq	s2limb
	subq	A_x, A_z
	leaq	(%rax), A_x
	movq	ADDR(s1p,index,eval(2*unroll_index+1)), %rax
	adcq	A_y, B_x
	movq	A_z, ADDR(rp,index,eval(2*unroll_index-2))
	movq	ADDR(rp,index,eval(2*unroll_index)), A_z
	adcq	$0, B_y
	leaq	(%rdx), A_y

	mulq	s2limb
	subq	B_x, B_z
	leaq	(%rax), B_x
	movq	ADDR(s1p,index,eval(2*unroll_index+2)), %rax
	adcq	B_y, A_x
	movq	B_z, ADDR(rp,index,eval(2*unroll_index-1))
	movq	ADDR(rp,index,eval(2*unroll_index+1)), B_z
	adcq	$0, A_y
	leaq	(%rdx), B_y
')

	addq	$UNROLL_SIZE, index
	jnz	L_mpn_submul_1_main_loop

L_mpn_submul_1_finish:
				C At this point we should have:
				C
				C index = n-1
				C
				C A_x = low_mul[index-2] + carry_in
				C A_y = high_mul[index-2] + CF
				C A_z = rp[index-2]
				C
				C B_x = low_mul[index-1]
				C B_y = high_mul[index-1]
				C B_z = rp[index-1]
				C
				C rax = s1p[index]
	mulq	s2limb
	subq	A_x, A_z
	movq	%rax, A_x
	movq	A_z, ADDR(rp,index,-2)
	movq	ADDR(rp,index,0), A_z
	adcq	A_y, B_x
	adcq	$0, B_y
	movq	%rdx, A_y
				C At this point we should have:
				C
				C index = n-1
				C
				C B_x = low_mul[index-1] + carry_in
				C B_y = high_mul[index-1] + CF
				C B_z = rp[index-1]
				C
				C A_x = low_mul[index]
				C A_y = high_mul[index]
				C A_z = rp[index]
	subq	B_x, B_z
	movq	B_z, ADDR(rp,index,-1)
	adcq	B_y, A_x
	adcq	$0, A_y
				C At this point we should have:
				C
				C index = n-1
				C
				C A_x = low_mul[index] + carry_in
				C A_y = high_mul[index] + CF
				C A_z = rp[index]
				C
	subq	A_x, A_z
	movq	A_z, ADDR(rp,index,0)
	adcq	$0, A_y

	movq	A_y, %rax
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
EPILOGUE()
