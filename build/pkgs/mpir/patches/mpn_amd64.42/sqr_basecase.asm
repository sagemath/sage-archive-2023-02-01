dnl  AMD K8 mpn_sqr_basecase -- square an mpn number.

dnl  This file is just an adaptation of similar file in the k7 directory.
dnl  Adapted by P. Gaudry in April 2005.
dnl  Here is the copyright of the original k7 version:

dnl  Copyright 1999, 2000, 2001, 2002 Free Software Foundation, Inc.
dnl
dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can rrdistribute it and/or
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

deflit(SQR_KARATSUBA_THRESHOLD_MAX, 34)

ifdef(`SQR_KARATSUBA_THRESHOLD_OVERRIDE',
 `define(`SQR_KARATSUBA_THRESHOLD',SQR_KARATSUBA_THRESHOLD_OVERRIDE)')

m4_config_gmp_mparam(`SQR_KARATSUBA_THRESHOLD')

C  deflit(UNROLL_COUNT, eval(SQR_KARATSUBA_THRESHOLD-3))

deflit(UNROLL_COUNT, 31)

C void mpn_sqr_basecase (mp_ptr dst, mp_srcptr src, mp_size_t size);
C
C With a SQR_KARATSUBA_THRESHOLD around 50 this code is about 1500 bytes,
C which is quite a bit, but is considered good value since squares big
C enough to use most of the code will be spending quite a few cycles in it.


define(param_dst, %rdi)
define(param_src, %rsi)
define(param_size, %rdx)

define(PARAM_DST, %r8)
define(PARAM_SRC, %r9)
define(PARAM_SIZE, %r10)

	TEXT
	ALIGN(32)
PROLOGUE(mpn_sqr_basecase)

	movq	param_size, %rcx
	movq	param_src, %rax
	cmpq	$2, %rcx

	movq	param_dst, %rdx
	je	L(two_limbs)
	ja	L(three_or_more)


C------------------------------------------------------------------------------
C one limb only
	C rax	src
	C rcx	size
	C rdx	dst
	C
	C rsi   src
	C rdi   dst

	movq	(%rsi), %rax
	mulq	%rax
	movq	%rdx, 8(%rdi)
	movq	%rax, (%rdi)
	ret


C------------------------------------------------------------------------------
C
C Using the read/modify/write "add"s seems to be faster than saving and
C restoring registers.  Perhaps the loads for the first set hide under the
C mul latency and the second gets store to load forwarding.

	ALIGN(16)
L(two_limbs):
	C rax
	C rbx
	C rcx
	C rdx
	C rsi	src
	C rdi   dst
	C
	C r8   s0
	C r9   s1

	movq	(%rsi), %r8
	movq	8(%rsi), %r9

        movq	%r8, %rax
	mulq	%rax		C src[0]^2

	movq	%rax, (%rdi)	C dst[0]
	movq	%r9, %rax

	movq	%rdx, 8(%rdi)	C dst[1]

	mulq	%rax		C src[1]^2

	movq	%rax, 16(%rdi)	C dst[2]
	movq	%r8, %rax

	movq	%rdx, 24(%rdi)	C dst[3]

	mulq	%r9		C src[0]*src[1]

 	addq	%rax, 8(%rdi)
	adcq	%rdx, 16(%rdi)
	adcq	$0, 24(%rdi)
	ASSERT(nc)

	addq	%rax, 8(%rdi)
	adcq	%rdx, 16(%rdi)
	adcq	$0, 24(%rdi)
	ASSERT(nc)

	ret


C------------------------------------------------------------------------------
defframe(SAVE_RBX,  -8)
defframe(SAVE_RBP, -16)
deflit(STACK_SPACE, 16)

L(three_or_more):
	subq	$STACK_SPACE, %rsp
	cmpq	$4, %rcx
	jae	L(four_or_more)
deflit(`FRAME',STACK_SPACE)


C------------------------------------------------------------------------------
C Three limbs
C
C Writing out the loads and stores separately at the end of this code comes
C out about 10 cycles faster than using adcls to memory.

	C rax	src
	C rcx	size
	C rdx	dst

	movq	%rbx, SAVE_RBX
	movq	%rax, %rbx	C src
	movq	(%rax), %rax

	movq	%rdx, %rcx	C dst

	mulq	%rax		C src[0] ^ 2

	movq	%rax, (%rcx)
	movq	8(%rbx), %rax
	movq	%rdx, 8(%rcx)

	mulq	%rax		C src[1] ^ 2

	movq	%rax, 16(%rcx)
	movq	16(%rbx), %rax
	movq	%rdx, 24(%rcx)

	mulq	%rax		C src[2] ^ 2

	movq	%rax, 32(%rcx)
	movq	(%rbx), %rax
	movq	%rdx, 40(%rcx)

	mulq	8(%rbx)		C src[0] * src[1]

	movq	%rax, %rsi
	movq	(%rbx), %rax
	movq	%rdx, %rdi

	mulq	16(%rbx)		C src[0] * src[2]

	addq	%rax, %rdi
	movq	%rbp, SAVE_RBP
	movq	$0, %rbp

	movq	8(%rbx), %rax
	adcq	%rdx, %rbp

	mulq	16(%rbx)		C src[1] * src[2]

	xorq	%rbx, %rbx
	addq	%rax, %rbp

	adcq	$0, %rdx

	C rax
	C rbx	zero, will be dst[5]
	C rcx	dst
	C rdx	dst[4]
	C rsi	dst[1]
	C rdi	dst[2]
	C rbp	dst[3]

	adcq	$0, %rdx
	addq	%rsi, %rsi

	adcq	%rdi, %rdi
	movq	8(%rcx), %rax

	adcq	%rbp, %rbp

	adcq	%rdx, %rdx

	adcq	$0, %rbx
	addq	%rax, %rsi
	movq	16(%rcx), %rax

	adcq	%rax, %rdi
	movq	24(%rcx), %rax
	movq	%rsi, 8(%rcx)

	adcq	%rax, %rbp
	movq	32(%rcx), %rax
	movq	%rdi, 16(%rcx)

	adcq	%rax, %rdx
	movq	40(%rcx), %rax
	movq	%rbp, 24(%rcx)

	adcq	%rbx, %rax
	ASSERT(nc)
	movq	SAVE_RBX, %rbx
	movq	SAVE_RBP, %rbp

	movq	%rdx, 32(%rcx)
	movq	%rax, 40(%rcx)
	addq	$FRAME, %rsp

	ret


C------------------------------------------------------------------------------
L(four_or_more):

C First multiply src[0]*src[1..size-1] and store at dst[1..size].
C Further products are added in rather than stored.

	C rax	src
	C rbx
	C rcx	size
	C rdx	dst
	C rsi
	C rdi
	C rbp

defframe(`VAR_COUNTER',-24)
defframe(`VAR_JMP',    -32)
deflit(EXTRA_STACK_SPACE, 16)
	movq	param_dst, PARAM_DST
	movq	param_src, PARAM_SRC
	movq	%rcx, PARAM_SIZE

	movq	%rbx, SAVE_RBX
	leaq	(%rdx,%rcx,8), %rdi	C &dst[size]

	movq	%rbp, SAVE_RBP
	leaq	(%rax,%rcx,8), %rsi	C &src[size]

	movq	(%rax), %rbp		C multiplier
	movq	$0, %rbx
	decq	%rcx

	negq	%rcx
	subq	$EXTRA_STACK_SPACE, %rsp
deflit(`FRAME', STACK_SPACE+EXTRA_STACK_SPACE)

L(mul_1):
	C rax	scratch
	C rbx	carry
	C rcx	counter
	C rdx	scratch
	C rsi	&src[size]
	C rdi	&dst[size]
	C rbp	multiplier

        movq    (%rsi,%rcx,8), %rax

        mulq    %rbp

        addq    %rbx, %rax
        movq    %rax, (%rdi,%rcx,8)
        movq    $0, %rbx

        adcq    %rdx, %rbx
        incq    %rcx
        jnz     L(mul_1)


C Add products src[n]*src[n+1..size-1] at dst[2*n-1...], for each n=1..size-2.
C
C The last two products, which are the bottom right corner of the product
C triangle, are left to the end.  These are src[size-3]*src[size-2,size-1]
C and src[size-2]*src[size-1].  If size is 4 then it's only these corner
C cases that need to be done.
C
C The unrolled code is the same as in mpn_addmul_1, see that routine for
C some comments.
C
C VAR_COUNTER is the outer loop, running from -size+4 to -1, inclusive.
C
C VAR_JMP is the computed jump into the unrolled code, stepped by one code
C chunk each outer loop.
C
C K7 does branch prrdiction on indirect jumps, which is bad since it's a
C different target each time.  There seems no way to avoid this.

dnl  This value also hard coded in some shifts and adds
deflit(CODE_BYTES_PER_LIMB, 25)

dnl  With the unmodified &src[size] and &dst[size] pointers, the
dnl  displacements in the unrolled code fit in a byte for UNROLL_COUNT
dnl  values up to 31, but above that an offset must be added to them.

deflit(OFFSET,
ifelse(eval(UNROLL_COUNT>15),1,
eval((UNROLL_COUNT-15)*8),
0))

dnl  Because the last chunk of code is generated differently, a label placed
dnl  at the end doesn't work.  Instead calculate the implied end using the
dnl  start and how many chunks of code there are.

deflit(UNROLL_INNER_END,
`L(unroll_inner_start)+eval(UNROLL_COUNT*CODE_BYTES_PER_LIMB)')

	C rax
	C rbx	carry
	C rcx
	C rdx
	C rsi	&src[size]
	C rdi	&dst[size]
	C rbp

	movq	PARAM_SIZE, %rcx
	movq	%rbx, (%rdi)

	subq	$4, %rcx
	jz	L(corner)

	negq	%rcx
ifelse(OFFSET,0,,`subq	$OFFSET, %rdi')
ifelse(OFFSET,0,,`subq	$OFFSET, %rsi')

	movq	%rcx, %rdx
	shlq	$4, %rcx
	movq	%rdx, %r11
	shlq	$3, %r11
	addq	%rdx, %r11
	C   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC Change here!!!!

ifdef(`PIC',`
	call	L(pic_calc)
L(here):
',`
	leaq	UNROLL_INNER_END-eval(2*CODE_BYTES_PER_LIMB)(%rcx,%r11), %rcx
')

	C The calculated jump mustn't come out to before the start of the
	C code available.  This is the limit UNROLL_COUNT puts on the src
	C operand size, but checked here directly using the jump address.
	ASSERT(ae,
	`movq_text_address(L(unroll_inner_start), %rax)
	cmpq	%rax, %rcx')


C------------------------------------------------------------------------------
	ALIGN(16)
L(unroll_outer_top):
	C rax
	C rbx	high limb to store
	C rcx	VAR_JMP
	C rdx	VAR_COUNTER, limbs, negative
	C rsi	&src[size], constant
	C rdi	dst ptr, high of last addmul
	C rbp

	movq	-24+OFFSET(%rsi,%rdx,8), %rbp	C next multiplier
	movq	-16+OFFSET(%rsi,%rdx,8), %rax	C first of multiplicand

	movq	%rdx, VAR_COUNTER

	mulq	%rbp

define(cmovX,`ifelse(eval(UNROLL_COUNT%2),0,`cmovz $@',`cmovnz $@')')

	testb	$1, %cl
C	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccc
	movq	%rdx, %rbx	C high carry
	movq	%rcx, %rdx	C jump

	movq	%rax, %rcx	C low carry
	cmovX(	%rbx, %rcx)	C high carry reverse
	cmovX(	%rax, %rbx)	C low carry reverse

	leaq	CODE_BYTES_PER_LIMB(%rdx), %rax
	xorq	%rdx, %rdx
	leaq	8(%rdi), %rdi

	movq	%rax, VAR_JMP

	jmp	*%rax


ifdef(`PIC',`
L(pic_calc):
	addq	(%rsp), %rcx
	addq	$UNROLL_INNER_END-eval(2*CODE_BYTES_PER_LIMB)-L(here), %rcx
	addq	%r11, %rcx
	ret
')


	C Must be an even address to preserve the significance of the low
	C bit of the jump address indicating which way around rcx/rbx should
	C start.
	ALIGN(2)

L(unroll_inner_start):
	C rax	next limb
	C rbx	carry high
	C rcx	carry low
	C rdx	scratch
	C rsi	src
	C rdi	dst
	C rbp	multiplier

forloop(`i', UNROLL_COUNT, 1, `
	deflit(`disp_src', eval(-i*8 + OFFSET))
	deflit(`disp_dst', eval(disp_src - 8))

	m4_assert(`disp_src>=-128 && disp_src<128')
	m4_assert(`disp_dst>=-128 && disp_dst<128')

ifelse(eval(i%2),0,`
.byte 0x90
        adcq    %rdx, %rbx
Zdisp(	movq,	disp_src,(%rsi), %rax)

        mulq	%rbp

Zdisp(  addq,	%rcx, disp_dst,(%rdi))
	movq	$0, %rcx

	adcq	%rax, %rbx

',`
	dnl  this bit comes out last
.byte 0x90
	adcq	%rdx, %rcx
Zdisp(  movq,	disp_src,(%rsi), %rax)

	mulq    %rbp

Zdisp(	addq,	%rbx, disp_dst,(%rdi))

ifelse(forloop_last,0,
`	movq	$0, %rbx')

	adcq    %rax, %rcx
')
')

	C rax	next limb
	C rbx	carry high
	C rcx	carry low
	C rdx	scratch
	C rsi	src
	C rdi	dst
	C rbp	multiplier

        adcq    $0, %rdx
	addq	%rcx, -8+OFFSET(%rdi)
	movq	VAR_JMP, %rcx

        adcq    $0, %rdx

	movq	%rdx, m4_empty_if_zero(OFFSET) (%rdi)
	movq	VAR_COUNTER, %rdx

	incq	%rdx
	jnz	L(unroll_outer_top)


ifelse(OFFSET,0,,`
	addq	$OFFSET, %rsi
	addq	$OFFSET, %rdi
')


C------------------------------------------------------------------------------
L(corner):
	C rsi	&src[size]
	C rdi	&dst[2*size-5]

	movq	-24(%rsi), %rbp
	movq	-16(%rsi), %rax
	movq	%rax, %rcx

	mulq	%rbp

	addq	%rax, -8(%rdi)
	movq	-8(%rsi), %rax

	adcq	$0, %rdx
	movq	%rdx, %rbx
	movq	%rax, %rsi

	mulq	%rbp

	addq	%rbx, %rax

	adcq	$0, %rdx
	addq	%rax, (%rdi)
	movq	%rsi, %rax

	adcq	$0, %rdx
	movq	%rdx, %rbx

	mulq	%rcx

	addq	%rbx, %rax
	movq	%rax, 8(%rdi)

	adcq	$0, %rdx
	movq	%rdx, 16(%rdi)



C Left shift of dst[1..2*size-2], high bit shifted out becomes dst[2*size-1].

L(lshift_start):
	movq	PARAM_SIZE, %rax
	movq	PARAM_DST, %rdi

	movq	%rax, %r11
	shlq	$1, %r11
	leaq	(%rdi,%r11,8), %rdi
	notq	%rax			C -size-1, preserve carry

	leaq	2(%rax), %rax		C -(size-1)
	movq	%rax, %r11
	shlq	$1, %r11

	xorq	%rcx, %rcx		C clear carry

L(lshift):
	C rax	counter, negative
	C rbx
	C rcx
	C rdx
	C rsi
	C rdi	dst, pointing just after last limb
	C rbp

	rclq	-8(%rdi,%r11,8)
	rclq	(%rdi,%r11,8)
	incq	%r11
	incq	%r11
	incq	%rax
	jnz	L(lshift)

	setc	%al

	movq	PARAM_SRC, %rsi
	movq	%rax, -8(%rdi)		C dst most significant limb

	movq	PARAM_SIZE, %rcx


C Now add in the squares on the diagonal, src[0]^2, src[1]^2, ...,
C src[size-1]^2.  dst[0] hasn't yet been set at all yet, and just gets the
C low limb of src[0]^2.

	movq	(%rsi), %rax		C src[0]

	mulq	%rax

	leaq	(%rsi,%rcx,8), %rsi	C src point just after last limb
	negq	%rcx

	movq 	%rcx, %r11
	shlq	$1, %r11

	movq	%rax, (%rdi,%r11,8)	C dst[0]
	incq	%rcx
	incq	%r11
	incq	%r11

L(diag):
	C rax	scratch
	C rbx	scratch
	C rcx	counter, negative
	C rdx	carry
	C rsi	src just after last limb
	C rdi	dst just after last limb
	C rbp

	movq	%rdx, %rbx
	movq	(%rsi,%rcx,8), %rax

	mulq	%rax

	addq	%rbx, -8(%rdi,%r11,8)
	adcq	%rax, (%rdi,%r11,8)
	adcq	$0, %rdx

	incq	%rcx
	incq    %r11
	incq    %r11

	jnz	L(diag)

	movq	SAVE_RBX, %rbx

	addq	%rdx, -8(%rdi)		C dst most significant limb

	movq	SAVE_RBP, %rbp
	addq	$FRAME, %rsp

	ret

EPILOGUE()
